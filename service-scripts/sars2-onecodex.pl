
=head1 NAME

    sars2-onecodex - OneCodex protocol for SARS2 assembly 
    
=head1 SYNOPSIS

    sars2-onecodex PE-read1 PE-read2 output-base output-dir
    sars2-onecodex SE-read output-base output-dir

=head1 DESCRIPTION

Assemble the input reads using the OneCodex protocol
L<https://github.com/onecodex/sars-cov-2>.    

Adapted from https://github.com/onecodex/sars-cov-2/blob/master/covid19_call_variants.sh

The generated FASTA consensus sequence is written to output-dir/output-base.fasta.

=cut

use strict;
use Getopt::Long::Descriptive;
use IPC::Run 'run';
use Bio::P3::SARS2Assembly qw(artic_reference artic_bed run_cmds reference_gff_path);
use Data::Dumper;
use File::Basename;
use File::Temp;

$ENV{PATH} = "$ENV{KB_RUNTIME}/samtools-1.11/bin:$ENV{KB_RUNTIME}/bcftools-1.9/bin:$ENV{PATH}";

my($opt, $usage) = describe_options("%c %o read1 [read2] output-base output-dir",
				    ["output-name|n=s" => "Output name for sequence (in the fasta file). Defaults to output-base"],
				    ["threads|j=i" => "Number of threads to use", { default => 1 }],
				    ["min-quality|q=i" => "Minimum read quality", { default => 20 }],
				    ["artic-version|a=i" => "ARTIC primer version", { default => 3}],
				    ["length-threshold|l=i" => "Max length to be considered short read sequencing", { default => 600 }],
				    ["keep-intermediates|k" => "Save all intermediate files"],
				    ["help|h"      => "Show this help message"],
				    );

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 3 || @ARGV == 4;

my $mode;

my($se_read, $pe_read_1, $pe_read_2);

my %artic_versions = (1 => 1, 2 => 1, 3 => 1);
if (!$artic_versions{$opt->artic_version})
{
    die "Invalid artic version\n";
}

my @inputs;
if (@ARGV == 3)
{
    $mode = "SE";
    $se_read  = shift;
    @inputs = ($se_read);
}
elsif (@ARGV == 4)
{
    $mode = "PE";
    $pe_read_1 = shift;
    $pe_read_2 = shift;
    @inputs = ($pe_read_1, $pe_read_2);
}
else
{
    die $usage->text;
}


#
# Set output directories
#

my $base = shift;
my $out_dir = shift;

my $int_dir;			# intermediate files

-d $out_dir or mkdir($out_dir) or die "Cannot create $out_dir: $!";

if ($opt->keep_intermediates)
{
    $int_dir = $out_dir;
}
else
{
    $int_dir = File::Temp->newdir(CLEANUP => 1);
}

$base =~ m,/, and die "Output base may not have slash characters\n";

my $output_name = $opt->output_name || $base;
$output_name =~ s/\s+/_/g;

#
# Mapping
# 

# mapping_mode=$(
#     head -n 400 "${input_fastq}" \
#     | awk \
#     -v "thresh=${length_threshold}" \
#     'NR % 4 == 2 {s+= length}END {if (s/(NR/4) > thresh) {print "map-ont"} else {print "sr"}}'
#    )

#
# Determine mapping mode (short or long) based on average read size of the start of the input
#

open(F, "-|", "gunzip", "-c", "-q",  "-d", "-f", $inputs[0]) or die "Cannot open $inputs[0]: $!";

my($count, $total);
while (<F>)
{
    last if $. > 400;
    if ($. % 4 == 2)
    {
	chomp;
	$count++;
	$total += length($_);
    }
}
close(F);
my $avg = $total / $count;
my $mapping_mode = ($avg > $opt->length_threshold) ? "map-ont" : "sr";

my @minimap_opts = (-K => "20M",
		    '-a',
		    -x => $mapping_mode,
		    -t => $opt->threads);

# Trim polyA tail for alignment (33 bases)
my $reference = artic_reference($opt->artic_version);
-f $reference or die "Cannot read reference $reference\n";
my $trimmed = "$int_dir/reference_trimmed.fa";

my $ok = run_cmds(["seqtk", "trimfq", "-e", 33, $reference], '>',  $trimmed);

$ok or die "Failure $? running seqtk\n";

#
# Run the mapper
# 

# shellcheck disable=SC2086
run_cmds(["minimap2", @minimap_opts,
	  $trimmed,
	  @inputs],
	 '|',
	 ["samtools",
	  "view",
	  "-u",
	  "-h",
	  "-q", $opt->min_quality,
	  "-F", 4,
	  "-"],
	 '|',
	 ["samtools",
	  "sort",
	  "--threads", $opt->threads,
	  "-"],
	 '>',
	 "$int_dir/$base.sorted.bam");

run_cmds(["samtools", "index", "$int_dir/$base.sorted.bam"]);

my $ivar_file = "$int_dir/$base.ivar";

run_cmds(["ivar",
	  "trim",
	  "-e",
	  "-q", 0,
	  "-i", "$int_dir/$base.sorted.bam",
	  "-b", artic_bed($opt->artic_version),
	  "-p", $ivar_file]);

run_cmds(["samtools",
	  "sort",
	  "$ivar_file.bam",
	  "--threads", $opt->threads],
	 ">",
	 "$int_dir/$base.isorted.bam");

run_cmds(["samtools", "index", "$int_dir/$base.isorted.bam"]);

run_cmds(["samtools",
	  "mpileup",
	  "--fasta-ref", $reference,
	  "--max-depth", 0,
	  "--count-orphans",
	  "--no-BAQ",
	  "--min-BQ", 0,
	  "$int_dir/${base}.isorted.bam"],
	 '>',
	 "$out_dir/$base.pileup");

run_cmds(["ivar",
	  "variants",
	  "-p", $ivar_file,
	  "-r", $reference,
	  "-g", reference_gff_path,
	  "-t", 0.6],
	 "<",
	 "$out_dir/$base.pileup");

run_cmds(["ivar",
	  "consensus",
	  "-p", $ivar_file,
	  "-m", 1,
	  "-t", 0.6,
	  "-n", "N"],
	 "<",
	 "$out_dir/$base.pileup");

run_cmds(["sed",
	  '/>/ s/$/ | One Codex consensus sequence/'],
	 "<",
	 "$ivar_file.fa",
	 '|',
	 [qw(seqtk seq -l 60 -)],
	 ">",
	 "$out_dir/$base.fasta");

run_cmds(["cat", $reference, "$out_dir/$base.fasta"],
	 '|',
	 ["mafft", "--auto", "-"],
	 '>',
	 "$out_dir/$base.align");

system("mv", "$int_dir/$base.isorted.bam", "$out_dir/$base.sorted.bam");
system("mv", "$int_dir/$base.isorted.bam.bai", "$out_dir/$base.sorted.bam.bai");
system("gzip", "$out_dir/$base.pileup");

run_cmds(["sars2-compute-ivar-vaf", "$ivar_file.tsv"],
	 '>',
	 "$out_dir/$base.variants.tsv");

__END__

my $samfile = "$int_dir/$base.sam";
my $bamfile = "$out_dir/$base.bam";
my $vcf = "$out_dir/$base.vcf";
my $vcf2 = "$out_dir/variants.vcf";
my $consensusfasta = "$out_dir/$base.fasta";

my $threads = $opt->threads;

#
# Adjust path to put the right version of tools in place.
#
$ENV{PATH} = "$ENV{KB_RUNTIME}/samtools-1.9/bin:$ENV{KB_RUNTIME}/bcftools-1.9/bin:$ENV{PATH}";

my $reference = "$int_dir/reference.fasta";
open(RIN, "<", $reference_base) or die "Cannot open reference $reference_base: $!";
open(ROUT, ">", $reference) or die "Cannot open reference $reference: $!";
{
    my $seen;
    while (<RIN>)
    {
	if (/^>/)
	{
	    if ($seen)
	    {
		die "unexpected multiple-contig reference found in $reference_base";
	    }
	    $seen = 1;
	    print ROUT ">$output_name\n";
	}
	else
	{
	    print ROUT $_;
	}
    }
    close(RIN);
    close(ROUT);
    run_cmds(["bowtie2-build", $reference, $reference]);
}

#
# Paired end and single end differ at the start. Single end mode is adapted from the
# CDC protocol.
#

my @min_rl = ("-m", $opt->minimum_read_length) if $opt->minimum_read_length;

if ($mode eq 'PE')
{
    #
    # Step 1. Adapter trimming.
    #
    
    my $trim1 = "$int_dir/$base.trim.1.fastq";
    my $trim2 = "$int_dir/$base.trim.2.fastq";
    
    my $t1 = 1;
    my $t2 = 1;
    if ($threads > 1)
    {
	$t1 = int($threads / 2);
	$t2 = $t2 = $threads - $t1;
    }
    
    my @cutadapt1 = qw(cutadapt
		       -g GTTTCCCAGTCACGATA
		       -G GTTTCCCAGTCACGATA
		       -a TATCGTGACTGGGAAAC
		       -A TATCGTGACTGGGAAAC
		       -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT 
		       -G ACACTCTTTCCCTACACGACGCTCTTCCGATCT
		       -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
		       -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
		       -n 3
		       -q 25);

    push(@cutadapt1,
	 "-j", $t1,
	 "--interleaved", $pe_read_1, $pe_read_2,
	 @min_rl,
	);
    
    my @cutadapt2 = qw(cutadapt
		       --interleaved
		       -u 30);
    push(@cutadapt2, 
	 "-j", $t2,
	 "-o", $trim1,
	 "-p", $trim2,
	 @min_rl,
	 "-");
    
    run_cmds(\@cutadapt1, '|', \@cutadapt2);
    
    #
    # Step 2. Mapping with bowtie. 
    #
    
    my @bowtie = ("bowtie2",
		  "--sensitive-local",
		  "-p", $threads,
		  "-x", $reference,
		  "-1", $trim1,
		  "-2", $trim2,
		  "-S", $samfile);
    run_cmds(\@bowtie);
}
elsif ($mode eq 'SE')
{
    #
    # Step 1. Adapter trimming.
    #
    
    my $trim = "$int_dir/$base.trim.fastq";
    
    my @cutadapt = qw(cutadapt
		  -g GTTTCCCAGTCACGATA
		  -a TATCGTGACTGGGAAAC
		  -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT 
		  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
		  -n 3
		  -q 25);
    push(@cutadapt,
	 "-j", $threads,
	 $se_read, "-o", $trim,
	 @min_rl,
	);

    run_cmds(\@cutadapt);
    
    #
    # Step 2. Mapping with bowtie. 
    #
    
    my @bowtie = ("bowtie2",
	      "--sensitive-local",
	      "-p", $threads,
	      "-x", $reference,
	      "-U", $trim,
	      "-S", $samfile);
    run_cmds(\@bowtie);
}
   
#
# Create bam format output
#

run_cmds(["samtools", "view", "-b", $samfile],
	 "|",
	 ["samtools", "sort", "-", "--threads", $threads, "-o", $bamfile]);

run_cmds(["samtools", "index", $bamfile]);

#
# Create consensus
#
my @con1 = ( qw(samtools mpileup -aa -d 8000 -uf), $reference, $bamfile);
my @con2 = qw(bcftools call -Mc);
my @con3 = ("tee", "-a", $vcf );
my @con4 = qw(vcfutils.pl vcf2fq -D 100000000);
push(@con4, "-d", $opt->min_depth);
my @con5 = qw(seqtk seq -A -);
my @con6 = qw(sed 2~2s/[actg]/N/g);
my @con7 = qw(seqtk seq -l 60 -);

run_cmds(\@con1, '|',
	 \@con2, '|',
	 \@con3, '|',
	 \@con4, '|',
	 \@con5, '|',
	 \@con6, '|',
	 \@con7, '>', $consensusfasta);

#
# Make sure any vcf files in the output folder are bgzipped/indexed.
#

for my $v (<$out_dir/*.vcf>)
{
    run_cmds(["bgzip", $v]);
}

for my $v (<$out_dir/*.vcf.gz>)
{
    if (! -f "$v.tbi")
    {
	run_cmds(["tabix", $v]);
    }
}

my $filter_params = "N_ALT == 0";
my @con1 = ( "bcftools", "filter", "-e", $filter_params, "$vcf.gz");
my @con2 = ("bgzip","-c");
run_cmds(\@con1, '|',
	\@con2, '>', "$vcf2.gz");

#
# Create coverage plot
#

run_cmds(["samtools", "depth", $bamfile], '>', "$out_dir/$base.depth");

if (-s "$out_dir/$base.depth")
{
    eval {
	$ENV{GDFONTPATH} = "/usr/share/fonts/liberation";
	my $plot = <<END;
set term png font "LiberationSans-Regular"
set xlabel "Position"
set ylabel "Depth"
set title "Coverage depth for $base"
set output "$out_dir/$base.png"
plot "$out_dir/$base.depth" using 2:3 with impulses title ""
set output
END
    run_cmds(["gnuplot"], "<", \$plot);
    };

    eval {
	$ENV{GDFONTPATH} = "/usr/share/fonts/liberation";
	my $plot = <<END;
set term png font "LiberationSans-Regular"
set yrange [0:250]
set xlabel "Position"
set ylabel "Depth"
set title "Coverage depth for $base"
set output "$out_dir/$base.detail.png"
plot "$out_dir/$base.depth" using 2:3 with impulses title ""
set output
END
    run_cmds(["gnuplot"], "<", \$plot);
    };


}

