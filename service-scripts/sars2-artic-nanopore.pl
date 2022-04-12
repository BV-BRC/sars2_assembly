
=head1 NAME

    sars2-artic-nanopore - Artic protocol for SARS2 assembly of Nanopore paired-end data
    
=head1 SYNOPSIS

    sars2-artic-nanopore reads output-base output-dir

=head1 DESCRIPTION

Assemble the input reads using the Artic protocol.

The generated FASTA consensus sequence is written to output-dir/output-base.fasta.

=cut

use strict;
use Getopt::Long::Descriptive;
use Bio::P3::SARS2Assembly 'run_cmds';
use File::Basename;
use File::Copy qw(move);
use File::Temp;

my($opt, $usage) = describe_options("%c %o reads output-base output-dir",
				    ['se-read|U=s' => "Single-end read file"],
				    ["output-name|n=s" => "Output name for sequence (in the fasta file). Defaults to output-base"],
				    ["scheme=s" => "Artic scheme", { default => "nCoV-2019/V3" }],
				    ["threads|j=i" => "Number of threads to use", { default => 1 }],
				    ["keep-intermediates|k" => "Save all intermediate files"],
				    ["help|h"      => "Show this help message"],
				    );

print($usage->text), exit 0 if $opt->help;

my $fastqfile;

#
# compatibility fix for new --se-read option.
#
if ($opt->se_read)
{
    $fastqfile = $opt->se_read;
    die($usage->text) unless @ARGV == 2;
}
else
{
    die($usage->text) unless @ARGV == 3;
    $fastqfile = shift;
}

my $base = shift;
my $out_dir = shift;

$base =~ m,/, and die "Output base may not have slash characters\n";

my $output_name = $opt->output_name || $base;
$output_name =~ s/\s+/_/g;

my $scheme_dir = Bio::P3::SARS2Assembly::artic_primer_schemes_path();

-d $out_dir or mkdir($out_dir) or die "Cannot create $out_dir: $!";
my $tmpdir = File::Temp->newdir(CLEANUP => 1);

my $threads = $opt->threads;

#
# Adjust path to put the right version of tools in place.
#
$ENV{PATH} = "$ENV{KB_RUNTIME}/samtools-1.9/bin:$ENV{KB_RUNTIME}/bcftools-1.9/bin:$ENV{PATH}";
$ENV{PATH} = "$ENV{KB_RUNTIME}/artic-ncov2019/bin:$ENV{PATH}";

#
# If not running with keep-intermediates, We run writing to our temp space,
# then move any files we are keeping into the destination.
#

my $run_dir = $opt->keep_intermediates ? $out_dir : $tmpdir;

my @artic = ("artic", "minion",
	     "--medaka",
	     "--normalise", 200,
	     "--threads", $threads,
	     "--scheme-directory", $scheme_dir,
	     "--read-file",  $fastqfile,
	     $opt->scheme,
	     "$run_dir/$base");
run_cmds(\@artic);

#
# We need to rewrite the output name here.
#

open(RIN, "<", "$run_dir/$base.consensus.fasta") or die "Cannot open consensus $run_dir/$base.consensus.fasta: $!";
open(ROUT, ">", "$run_dir/$base.fasta") or die "Cannot open output consensus $run_dir/$base.fasta: $!";
{
    my $seen;
    while (<RIN>)
    {
	if (/^>/)
	{
	    if ($seen)
	    {
		die "unexpected multiple-contig reference found in $base.consensus.fasta";
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
    unlink("$run_dir/$base.consensus.fasta");
}

if (!$opt->keep_intermediates)
{
    #
    # We keep the generated fasta, bam, and vcf files.
    #
    opendir(D, $tmpdir) or die "Cannot opendir $tmpdir: $!";
    while (my $f = readdir(D))
    {
	next if $f =~ /^$base.*\.nCoV-2019/;
	if ($f eq "$base.fasta" ||
	    $f =~ /minion\.log|\.bam|\.vcf\.gz/)
	{
	    move("$tmpdir/$f", "$out_dir/$f");
	}
    }
}

