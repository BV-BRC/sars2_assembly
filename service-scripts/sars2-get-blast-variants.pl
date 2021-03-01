#!/usr/bin/env perl
use strict;
use IPC::Run 'run';
#use XML::Simple;
use JSON::XS;
use File::Slurp;
use Data::Dumper;
use Getopt::Long;
use gjoseqlib;

my $usage = 'get-blast-variants.pl -s subject_seqs.dna -q spike_protein.aa


        -b blast program  [chose tblastn, blastp, blastn, or blastx], required                
        -q query sequence  required
        -s subject sequence for blasting against, required
        
        -l min hsp length (d = 50)
        -p min hsp percent identity d = (0.85)
        -t [file name] print table of blast match statistics
        -c [file name] print a consensus sequence of the recovered protein
        -g "gap character" d = "X", this is the character that gets inserted when
            multiple HSPS are spliced together. Blast fills in gaps in a match with 
            the X character.  
        -k keep the bastdb files [default = remove]
        -f strip fasta header after the first white space, default keeps whole string. 
        -r show frame [only works for tblastn]


        #these are currently turned off
        Blast Options
        #-threads = number of threads for the blast d = 8
        -dbsize  = -max_file_sz size limit of the makeblastdb d=10GB

        ---------
        Standard Output is:
        Subject ID\tFrames\tInsertions\tVariants
        Where
        -Frames are a comma separated list of all the frames that are found in the HSPs 
        -Insertions are defined as a the query position and the subject AAs at that postion e.g., 21-MMMMMMMM,45-G,68-C,91-M
        - Variants are one of the following:
            D614G  -- amino acid substitution where it is [query_position_subject]
            E281X  -- amino acid is an X-character in subject [usually happens when there are N nucleotides]
            A783@  -- where "@" is any symbol defined with the -g flag, these are missing regions in the query
            where the subject got spliced together.
            S112*  -- where the subject has a stop codon
            H245-  -- where the subject sequence has a credible gap         

';

my ($s_file, $q_file, $help, $table, $consensus, $keepdbs, $strip_header, $showframe, $program);
my $min_pid = 0.85; 
my $min_len = 50;
#my $threads = 8;
my $dbsize = "1GB"; 
my $gapc = "X";
my $opts = GetOptions( 'b=s'       => \$program,
		      'q=s'       => \$q_file,
		      's=s'       => \$s_file, 
		      'l=s'       => \$min_len,
		      'p=s'       => \$min_pid,
		      't=s'       => \$table,
		      'c=s'       => \$consensus,
		      'g=s'       => \$gapc,
		      'k'         => \$keepdbs,
		      'f'         => \$strip_header,
		      'r'         => \$showframe,
		      'dbsize=s'  => \$dbsize,
		      #'threads=i' => \$threads,
		      'h'         => \$help,);

if ($help){die "$usage\n";}
unless ($s_file ){die "must declare a subject file with -s \n\n$usage\n";}
unless ($q_file ){die "must declare a query file with -q \n\n$usage\n";}
unless($program){die "must declare a blast program with -b [tblastn or blastp]\n\n$usage\n";}

my $dbtype;
if ($program =~ /blastn/){$dbtype = "nucl";}
else{$dbtype = "prot";}


#Have to read the sequence to piece together multiple hsps
open (IN, "<$q_file"), or die "cannot open query file: $q_file\n"; 
my @full_query = &gjoseqlib::read_fasta(\*IN); 
my @full_query_seq = split ('', $full_query[0][2]);   
close IN; 

open (IN, "<$s_file"), or die "cannot open subj file: $s_file\n"; 
my @subseqs = &gjoseqlib::read_fasta(\*IN);  
my $max = scalar @subseqs;
print STDERR "\n BLASTing $max subject sequences\n\n"; 
close IN; 



my $ok = run(["makeblastdb", "-dbtype", $dbtype, "-in", $s_file, "-max_file_sz", $dbsize],
	     '>', '/dev/null');
   
$ok or die "Failure $? running makeblastdb";

my @cmd = ($program,
	   "-outfmt", "15",
	   "-subject", $s_file,
	   "-max_target_seqs", $max);


if (($program eq "blastn") || ($program eq "blastx"))
{
    push(@cmd, 
	 "-dust", "no",);
}
else
{
    push(@cmd,
	 '-seg', 'no');
}

my($out, $err);

print STDERR "Run: @cmd\n";
my $ok = run(\@cmd, '<', $q_file, '>', \$out, '2>', \$err);

$ok or die "Failure $? running @cmd\n";

if ($table)			#If you want an M8-like blast output
{
    open (TBL, ">$table");  
    print TBL "QID\tSID\tHSPs\tFrame\t%ID\tE\tBit\tQLen\tAli_len\tGaps\tQstart\tQend\tSstart\tSend\n";
}

if ($consensus)	 #If you want the recovered prot sequences of the subj
{	
    open (CON, ">$consensus");
}

my $blast = decode_json($out);

my $bl2seq  = $blast->{BlastOutput2}->[0]->{report}->{results}->{bl2seq};	
for my $i (0..$#{$bl2seq})
{
    my $qid     = $bl2seq->[$i]->{query_title};
    my $qlen    = $bl2seq->[$i]->{query_len};
    my $hits    = $bl2seq->[$i]->{hits};

    if (@$hits == 0)
    {
	print STDERR "No hits found\n";
	exit(2);
    }
    
    for my $j (0..$#{$hits})
    {
	#I do not understand the conditions where there are multiple subject
	#descriptions, if these arise, the zero below needs to be dealt with.
	my $sid   = $hits->[$j]->{description}->[0]->{title};  
	my $hsps  = $hits->[$j]->{hsps};	
	
	my $hash = {};
	my %dashes;
	my %frames;
	my @snps;
	
	for my $k (0..$#{$hsps})
	{
	    my $hsp_num  = $hsps->[$k]->{num};
	    my $ali_len  = $hsps->[$k]->{align_len};
	    my $id_len   = $hsps->[$k]->{identity};
	    my $q_start  = $hsps->[$k]->{query_from};
	    my $q_end    = $hsps->[$k]->{query_to};
	    my $s_start  = $hsps->[$k]->{hit_from};
	    my $s_end    = $hsps->[$k]->{hit_to};
	    my $mid_seq  = $hsps->[$k]->{midline};
	    my $q_seq    = $hsps->[$k]->{qseq};
	    my $s_seq    = $hsps->[$k]->{hseq};
	    my $eval     = $hsps->[$k]->{evalue};
	    my $gaps     = $hsps->[$k]->{gaps};
	    my $bit      = $hsps->[$k]->{bit_score};		
	    my $frame    = $hsps->[$k]->{hit_frame};	
	    my $positive = $hsps->[$k]->{positive};	
	    my $scr      = $hsps->[$k]->{score};	
	    
	    my $pid      = ($id_len/$ali_len);
	    my $inc      = $q_start;
	    my @query    = split ('', $q_seq);
	    my @subject  = split ('', $s_seq);
	    if($strip_header){ $sid =~ s/\s.+//g; $qid =~ s/\s.+//g; }
	    
	    #Filter HSPs on user input parameters:
	    if (($ali_len >= $min_len) && ($pid >= $min_pid))
	    {
		$frames{$frame} = 1;
		if ($table)
		{
		    print TBL "$qid\t$sid\t$hsp_num\t$frame\t$pid\t$eval\t$bit\t$qlen\t$ali_len\t$gaps\t$q_start\t$q_end\t$s_start\t$s_end\n";
		}		
		
		for my $l (0..$#query)
		{	
		    if ($query[$l] !~ /\-/)
		    {
			$hash->{$inc}->{'Q'} = $query[$l];
			$hash->{$inc}->{'S'} = $subject[$l];
			$inc ++; 
		    }
		    else     
		    {
			#populate the hash of insertions
			push @{$dashes{$inc}},$subject[$l];
		    }
		}
	    }
	}		
	
	#Merge HSPs, reconstruct the gaps in the subject protein sequence
	my $outlen = $qlen;
	if ($program eq "blastx"){$outlen = ($qlen/3);} # if you have a stop codon it will be represented as a character.
	
	for my $m (0..($outlen - 1))
	{
	    my $pos = ($m + 1);	#increment the position by one to match blast locs
	    if (exists $hash->{$pos})
	    {
		my $q_char = uc($hash->{$pos}->{Q});
		my $s_char = uc($hash->{$pos}->{S});
		if ($q_char ne $s_char)
		{
		    push @snps, "$q_char"."$pos"."$s_char";
		}
	    }
	    else
	    {
		my $seqchar = $full_query_seq[$m]; 
		push @snps, "$seqchar"."$pos"."$gapc"; #<-- this is a character that was missed by the blast
		# It is probably a gap in the subject, but there could be edge cases where it is not. 
				#populate the hash to build a consensus
		$hash->{$pos}->{Q} = $seqchar;
		$hash->{$pos}->{S} = $gapc;
	    }
	}
	
	# print stdout snp table
	my $nframes = keys %frames;
	print "$sid\t";
	if ($showframe) {print join (",", keys %frames), "\t";}
	
	my @ins;
	foreach (sort {$a <=> $b} keys %dashes)
	{
	    my $run = join ("", @{$dashes{$_}});; 
	    push @ins, "$_"."-"."$run"; 
	}
	print join (",", @ins), "\t"; 
	print join (",", @snps), "\n"; 
	
	#Build optional subject sequence for output
	if ($consensus)  
	{
	    my $seq_string;
	    foreach (sort {$a <=> $b} keys %{$hash})
	    {
		my $pos = $_;
		if (exists $dashes{$pos})
		{
		    my $ins = join ("", @{$dashes{$pos}});
		    $seq_string .= $ins
		    }
		$seq_string .= $hash->{$pos}->{S};
	    }
	    my $header = $sid;
	    &gjoseqlib::print_alignment_as_fasta(\*CON, [$header, "Query=$qid", $seq_string]);
	}
    }
}

#cleanup
unless ($keepdbs)
{
    if ($dbtype eq "nucl")
    {
	unlink ("$s_file.nin", "$s_file.nsq", "$s_file.nhr");
    }
    else
    {
	unlink ("$s_file.pin", "$s_file.phr", "$s_file.psq");
    }
}










