
=head1 NAME

    sars2-cdc-nanopore - CDC protocol for SARS2 assembly of Nanopore data
    
=head1 SYNOPSIS

    sars2-cdc-nanopore readfile output-base output-dir

=head1 DESCRIPTION

Assemble the input reads using the Nanopore assembly protocol
from the CDC "Protocols for SARS-C-V-2 sequencing" document.

The generated FASTA consensus sequence is written to output-dir/output-base.fasta.

=cut

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use Bio::P3::SARS2Assembly;
use File::Basename;
use File::Temp;
use Bio::P3::CmdRunner;
use JSON::XS;

my($opt, $usage) = describe_options("%c %o output-base output-dir",
				    ['se-read|U=s@' => "Single-end read file(s)", { default => [] }],
				    ["output-name|n=s" => "Output name for sequence (in the fasta file). Defaults to output-base"],
				    ["threads|j=i" => "Number of threads to use", { default => 1 }],
				    ["min-depth|d=i" => "Minimum depth required to call bases in consensus", { default => 3 }],
				    ["keep-intermediates|k" => "Save all intermediate files"],
				    ["help|h"      => "Show this help message"],
				    );

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 2;

my $runner = Bio::P3::CmdRunner->new;

my @reads = @{$opt->se_read};
if (@reads == 0)
{
    die "$0: No reads specified\n";
}

my $base = shift;
my $out_dir = shift;

$base =~ m,/, and die "Output base may not have slash characters\n";

my $output_name = $opt->output_name || $base;
$output_name =~ s/\s+/_/g;

my $reference_base = Bio::P3::SARS2Assembly::reference_fasta_path();
my $primer = Bio::P3::SARS2Assembly::primer_bedpe_path();

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
}
	

my $samfile = "$int_dir/$base.sam";
my $bamfile = "$out_dir/$base.bam";
my $vcf = "$out_dir/$base.vcf";
my $consensusfasta = "$out_dir/$base.fasta";
my $primerclippedbamfile = "$out_dir/$base.primerclipped.bam";
my $hdf = "$primerclippedbamfile.hdf";

my $threads = $opt->threads;

#
# Adjust path to put the right version of tools in place.
#
$ENV{PATH} = "$ENV{KB_RUNTIME}/samtools-1.9/bin:$ENV{KB_RUNTIME}/bcftools-1.9/bin:$ENV{PATH}";
$ENV{PATH} = "$ENV{KB_RUNTIME}/medaka/bin:$ENV{PATH}";

#
# Step 1. Adapter trimming.
#

my @cutadapt = qw(cutadapt
		  -m 300
		  -M 1200
		  -q 15 );

my @filtered;
for my $fastqfile (@reads)
{
    my $bn = basename($fastqfile);
    $bn =~ s/\.gz$//;
    $bn =~ s/\.[^.]+$//;
    my $fastqfiltered = "$int_dir/$bn.filtered.fastq";

    my @cmd = (@cutadapt, 
	       "-j", $threads,
	       "-o", $fastqfiltered,
	       $fastqfile);

    $runner->run(\@cmd);
    push(@filtered, $fastqfiltered);
}

#
# Step 2. Mapping with minimap2
#

my @minimap = ("minimap2",
	       "-L",
	       "-a",
	       "-x", "map-ont",
	       "-t", $threads,
	       $reference,
	       @filtered);
$runner->run(\@minimap, ">", $samfile);

#
# Create bam format output
#

$runner->run(["samtools", "view", "-b", $samfile],
	 "|",
	 ["samtools", "sort", "-", "--threads", $threads, "-o", $bamfile]);

$runner->run(["samtools", "index", $bamfile]);

#
# Primer clipping
#

my @clip = ("bamclipper.sh",
	    "-b", $bamfile,
	    "-p", $primer,
	    "-n", $threads,
	    "-o", $out_dir,
	    "-u", 80,
	    "-d", 80);
$runner->run(\@clip);

#
# Create consensus
#
my @con = ("medaka", "consensus",
	   "--model", "r941_min_high_g344",
	   "--threads", $threads,
	   $primerclippedbamfile,
	   $hdf);

$runner->run(\@con);

$runner->run(["medaka", "variant", $reference, $hdf, $vcf]);

my @mask = ('vcf_mask_lowcoverage.pl',
	    "--bam", $primerclippedbamfile,
	    "--reference", $reference,
	    "--vcf", $vcf,
	    "--refout", "$out_dir/$base.reference.masked.fasta",
	    "--consout", $consensusfasta,
	    "--depth", $opt->min_depth,
	    "--qual", 40);
$runner->run(\@mask);

#
# Make sure any vcf files in the output folder are bgzipped/indexed.
#

for my $v (<$out_dir/*.vcf>)
{
    $runner->run(["bgzip", $v]);
}

for my $v (<$out_dir/*.vcf.gz>)
{
    if (! -f "$v.tbi")
    {
	$runner->run(["tabix", $v]);
    }
}

#
# Create coverage plot
#

$runner->run(["samtools", "depth", $primerclippedbamfile], '>', "$out_dir/$base.depth");


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
    $runner->run(["gnuplot"], "<", \$plot);
    };
}

print STDERR  JSON::XS->new->pretty(1)->canonical(1)->encode($runner->report);
