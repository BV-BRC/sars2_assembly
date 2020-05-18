
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

my $threads = $opt->threads;

#
# Adjust path to put the right version of tools in place.
#
$ENV{PATH} = "$ENV{KB_RUNTIME}/artic/bin:$ENV{PATH}";

my @artic = ("artic", "minion",
	     "--medaka",
	     "--normalise", 200,
	     "--threads", $threads,
	     "--scheme-directory", $scheme_dir,
	     "--read-file",  $fastqfile,
	     $opt->scheme,
	     "$out_dir/$base");
run_cmds(\@artic);
	    
rename("$out_dir/$base.consensus.fasta", "$out_dir/$base.fasta");
