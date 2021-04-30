#
# Given a GTO, compute variation using
#
# Pangolin
# Jim's sars2-get-blast-variants
# (When ready) Maulik's alignment-based tool
#

use gjoseqlib;
use strict;
use Data::Dumper;
use GenomeTypeObject;
use Bio::P3::SARS2Assembly qw(reference_spike_aa_path mpath add_variants_to_gto add_quality_estimate_to_gto);
use Getopt::Long::Descriptive;
use IPC::Run qw(run);

my($opt, $usage) = describe_options("%c %o",
				    ["variants=s" => "Add variants.tsv from assembly"],
				    ["input|i=s" => "Input file"],
				    ["output|o=s" => "Output file"],
				    ["debug|d" => "Enable debugging"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 0;

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";

my $contigs = $genome_in->extract_contig_sequences_to_temp_file();

my $sequence_features = load_sequence_features(mpath . "/spike.sf");
my $loc = load_loc(mpath . "/spike.loc");

eval {
    if ($opt->variants)
    {
	add_variants_to_gto($opt->variants, $genome_in);
    }
    pangolin($genome_in, $contigs);
    get_blast_variants($genome_in, $contigs);
    maulik_mafft($genome_in, $sequence_features, $loc);
    add_quality_estimate_to_gto($genome_in);
};

if ($@)
{
    warn "Evaluation failure: $@";
}

unlink($contigs);

$genome_in->destroy_to_file($opt->output, { canonical => 1 });

sub pangolin
{
    my($gto, $contigs_file) = @_;

    my $out = File::Temp->new();
    close($out);
    my @cmd = ("pangolin", "--outfile", $out, $contigs_file);
    my $ok = run(\@cmd);
    if (!$ok)
    {
	warn "Pangolin run failed with $?\n";
	return;
    }

    if (!open(LIN, "<", $out))
    {
	warn "Cannot open pangolin lineage output $out: $!";
	return;
    }

    my $hdr = <LIN>;
    chomp $hdr;
    my @keys = split(/,/, $hdr);
    my %idx;
    for my $i (0..$#keys)
    {
	$idx{$keys[$i]} = $i;
    }
    my $data = <LIN>;
    chomp $data;
    my @data = split(/,/, $data);
    my $lin = $data[$idx{lineage}];
    my $prob = $data[$idx{probability}] if defined($idx{probability});
    my $conflict = $data[$idx{conflict}] if defined($idx{conflict});
    my $status = $data[$idx{status}];
    my $learn_vers = $data[$idx{pangoLEARN_version}];
    my $note = $data[$idx{note}];
    
    #
    # look up pangolin version
    #
    my $p_vers;
    my $ok = run(["pangolin", "-v"], ">", \$p_vers);
    if ($ok)
    {
	chomp $p_vers;
	$p_vers =~ s/^pang\S+\s+//;
    }

    print "lin=$lin prob=$prob status=$status\n";

    my $tool_md = {
	pangoLEARN_version => $learn_vers,
	pangolin_version => $p_vers,
    };

    my $var = {
	tool => 'pangolin',
	tool_metadata => $tool_md,
	variants => [],
	lineage => $lin,
	probability => $prob,
	conflict => $conflict,
	status => $status,
	notes => $note,
    };

    push(@{$gto->{computed_variants}}, $var);
}

sub get_blast_variants
{
    my($gto, $contigs_file) = @_;

    my $spike = reference_spike_aa_path;

    open(R, "<", reference_spike_aa_path) or die "Cannot open spike reference: $!";
    my($ref_id) = read_next_fasta_seq(\*R);
    close(R);

    my($out, $err);
    
    my @cmd = ("sars2-get-blast-variants",
	       "-b", "tblastn",
	       "-r",
	       "-q", $spike,
	       "-s", $contigs_file);
    print "@cmd\n";
    my $ok = run(\@cmd, ">", \$out, '2>', \$err);
    my $nohits;
    if (!$ok)
    {
	my $rc = ($? >> 8);
	if ($rc == 2)
	{
	    warn "No hits found\n";
	    $nohits = 1;
	}
	else
	{
	    warn "Error $? running @cmd\n$err\n";
	    return;
	}
    }

    my $vlist = [];
    my $tv = {
	tool => 'sars2-get-blast-variants',
	variants => $vlist,
    };

    open(M, "<", \$out);
    while (<M>)
    {
	print ;
	chomp;
	my($id, $frame, $ins, $var) = split(/\t/);
	my @ins = split(/,/, $ins);
	my @vars = split(/,/, $var);

	my $snps = [];
	for my $var (@vars)
	{
	    my($ref,$pos,$alt) = $var =~ /^(\D+)(\d+)(\D+)$/;
	    push(@$snps, {
		ref_aa => $1,
		feature_pos => int($2),
		alt_aa => $3,
	    });
	}
	#
	# NB this tool always looks at spike; gene is hardcoded here
	#
	push(@$vlist, {
	    reference => $ref_id,
	    gene => "S",
	    frame => $frame,
	    snps => $snps,
	});
    }
    push(@{$gto->{computed_variants}}, $tv);
    
}

sub maulik_mafft
{
    my($gto, $sfs, $locs)= @_;

    open(R, "<", reference_spike_aa_path) or die "Cannot open spike reference: $!";
    my($ref_id, $ref_def, $ref_seq) = read_next_fasta_seq(\*R);

    my @spike_prots;
    for my $feat ($gto->features)
    {
	push @spike_prots, $feat if $feat->{function} =~ /^(putative\s+)?surface\s+glycoprotein/;
    }

    if (@spike_prots != 1)
    {
	warn "Skipping evaluation due to number of spike proteins =" . scalar(@spike_prots) . "\n";
	return;
    }

    my $tmp = File::Temp->new();
    for my $feat (@spike_prots)
    {    
	print $tmp ">$feat->{id}\n$feat->{protein_translation}\n";
    }
    close($tmp);

    my $aln;
    my $err;
    my $ok = run(["mafft",
		  "--auto", "--keeplength",
		  "--addfragments", $tmp,
		  reference_spike_aa_path],
		 ">", \$aln,
		 '2>', \$err);

    open(ALN, "<", \$aln);
    my $out;

    my @refaa = split //, $ref_seq;

    my @out_vars;
    my $out_lin;
    
    my $seq;
    my $total;
    my @snps;
    while (my($seq_id, $seq_def, $seq) = read_next_fasta_seq(\*ALN))
    {
	if ($seq_id eq $ref_id)
	{
	    print STDERR "Skip ref $seq_id\n";
	    next;
	}
	$total++;

	my @aa = split //, uc $seq;
	my @vars = ();
	my @mysfs = ();
	for (my $i = 0; $i < length($seq)-1; $i++){
	    my $pos = $i+1;
	    if ($aa[$i] ne $refaa[$i] && $aa[$i] ne "X")
	    {
		push(@snps, {
		    ref_aa => $refaa[$i],
		    alt_aa => $aa[$i],
		    feature_pos => $pos,
		});
		push @vars, "$refaa[$i]$pos$aa[$i]";

		if ($sfs->{$pos})
		{
		    push @mysfs, split /,/, $sfs->{$pos};
		}
	    }
	}
	my %htmp;
	foreach my $sf (@mysfs){$htmp{$sf}=1;}
	@mysfs = keys %htmp;
	
	push(@out_vars, @vars);
	my $lineage = join(',', @vars);
	my $sequence_features = join(',', @mysfs);

	
	my $loc = "";
	foreach my $label (sort keys %$locs){
	    my $match = 1;
	    foreach my $var (split /,/, $locs->{$label})
	    {
		$match = 0 unless $lineage=~/$var/;
	    }
	    $loc = $label if $match == 1 && $loc eq "";
	}
	$out_lin = $loc;
	    
	    # print STDERR "$total\t$id\n";
    }
    close ALN;

    #
    # NB this tool always looks at spike; gene is hardcoded here
    #
    my $tv = {
	tool => 'variantTracker',
	variants => [ {
	    reference => $ref_id,
	    gene => "S",
	    snps => \@snps,
	} ],
	lineage => $out_lin,
    };
    push(@{$gto->{computed_variants}}, $tv);
    
}

sub load_sequence_features
{
    my($file) = @_;
    open(SF, "<", $file) or die "Cannot load sequence features from $file: $!";

    my $sfs = {};
    while (my $line=<SF>)
    {
	chomp $line;
	my ($start, $end, $label) = split /\t/, $line;
	for (my $pos=$start; $pos<=$end; $pos++)
	{
	    $sfs->{$pos} .= "$label,"; 
	} 
    }
    close(SF); 
    
    foreach my $pos (keys %$sfs){
	$sfs->{$pos}=~s/,$//;
    }
    return $sfs;
    
}

sub load_loc
{
    my($file) = @_;
    
    open(LOC, "<", $file) or die "Cannot read $file: $!";

    my $locs = {};
    
    while (my $line=<LOC>)
    {
	chomp $line;
	my ($label, $vars) = split /\t/, $line;
	$vars =~ s/\s*//g;
	$locs->{$label} = $vars;
    }
    close LOC; 
    return $locs;
}
