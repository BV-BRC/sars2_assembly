
=head1 NAME

    sars2-cdc-illumina - CDC protocol for SARS2 assembly of Illumina paired-end data
    
=head1 SYNOPSIS

    sars2-cdc-illumina read1 read2 output-base output-dir

=head1 DESCRIPTION

Assemble the input reads using the Illumina assembly protocol
from the CDC "Protocols for SARS-C-V-2 sequencing" document.

The generated FASTA consensus sequence is written to output-dir/output-base.fasta.

=cut

use strict;
use Getopt::Long::Descriptive;
use Bio::P3::SARS2Assembly 'run_cmds';
use File::Basename;
use File::Temp;

my($opt, $usage) = describe_options("%c %o read1 read2 output-base output-dir",
				    ["output-name|n=s" => "Output name for sequence (in the fasta file). Defaults to output-base"],
				    ["threads|j=i" => "Number of threads to use", { default => 1 }],
				    ["min-depth|d=i" => "Minimum depth required to call bases in consensus", { default => 100 }],
				    ["keep-intermediates|k" => "Save all intermediate files"],
				    ["help|h"      => "Show this help message"],
				    );

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 4;

my $read1 = shift;
my $read2 = shift;
my $base = shift;
my $out_dir = shift;

$base =~ m,/, and die "Output base may not have slash characters\n";

my $output_name = $opt->output_name || $base;
$output_name =~ s/\s+/_/g;

my $reference_base = Bio::P3::SARS2Assembly::reference_fasta_path();

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

my $trim1 = "$int_dir/$base.trim.1.fastq";
my $trim2 = "$int_dir/$base.trim.2.fastq";

my $samfile = "$int_dir/$base.sam";
my $bamfile = "$out_dir/$base.bam";
my $vcf = "$out_dir/$base.vcf";
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
# Step 1. Adapter trimming.
#

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
		  -m 75
		  -q 25);
push(@cutadapt1,
     "-j", $t1,
     "--interleaved", $read1, $read2
     );

my @cutadapt2 = qw(cutadapt
		   --interleaved
		   -m 75
		   -u 30);
push(@cutadapt2, 
     "-j", $t2,
     "-o", $trim1,
     "-p", $trim2,
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
}

