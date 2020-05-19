
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
use Getopt::Long::Descriptive;
use Bio::P3::SARS2Assembly 'run_cmds';
use File::Basename;
use File::Temp;

my($opt, $usage) = describe_options("%c %o reads output-base output-dir",
				    ["threads|j=i" => "Number of threads to use", { default => 1 }],
				    ["keep-intermediates|k" => "Save all intermediate files"],
				    ["help|h"      => "Show this help message"],
				    );

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 3;

my $fastqfile = shift;
my $base = shift;
my $out_dir = shift;

my $reference = Bio::P3::SARS2Assembly::reference_fasta_path();
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

my $fastqfiltered = "$int_dir/$base.filtered.fastq";

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
$ENV{PATH} = "$ENV{KB_RUNTIME}/medaka/bin:$ENV{PATH}";

#
# Step 1. Adapter trimming.
#

my @cutadapt = qw(cutadapt
		  -m 300
		  -M 1200
		  -q 15 );
push(@cutadapt, 
     "-j", $threads,
     "-o", $fastqfiltered,
     $fastqfile);

run_cmds(\@cutadapt);

#
# Step 2. Mapping with minimap2
#

my @minimap = ("minimap2",
	       "-L",
	       "-a",
	       "-x", "map-ont",
	       "-t", $threads,
	       $reference,
	       $fastqfiltered);
run_cmds(\@minimap, ">", $samfile);

#
# Create bam format output
#

run_cmds(["samtools", "view", "-b", $samfile],
	 "|",
	 ["samtools", "sort", "-", "--threads", $threads, "-o", $bamfile]);

run_cmds(["samtools", "index", $bamfile]);

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
run_cmds(\@clip);

#
# Create consensus
#
my @con = ("medaka", "consensus",
	   "--model", "r941_min_high_g344",
	   "--threads", $threads,
	   $primerclippedbamfile,
	   $hdf);

run_cmds(\@con);

run_cmds(["medaka", "variant", $reference, $hdf, $vcf]);

my @mask = ('vcf_mask_lowcoverage.pl',
	    "--bam", $primerclippedbamfile,
	    "--reference", $reference,
	    "--vcf", $vcf,
	    "--refout", "$out_dir/$base.reference.masked.fasta",
	    "--consout", $consensusfasta,
	    "--depth", 20,
	    "--qual", 40);
run_cmds(\@mask);

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

