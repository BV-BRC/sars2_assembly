
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
				    ["scheme=s" => "Artic scheme", { default => "nCoV-2019/V3" }],
				    ["threads|j=i" => "Number of threads to use", { default => 1 }],
				    ["keep-intermediates|k" => "Save all intermediate files"],
				    ["help|h"      => "Show this help message"],
				    );

print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 3;

my $fastqfile = shift;
my $base = shift;
my $out_dir = shift;

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

rename("$run_dir/$base.consensus.fasta", "$run_dir/$base.fasta");

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

