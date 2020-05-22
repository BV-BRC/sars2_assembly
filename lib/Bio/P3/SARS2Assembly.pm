package Bio::P3::SARS2Assembly;

#
# Mostly empty package; used as a placeholder 
#

use strict;
use Module::Metadata;
use IPC::Run 'run';
use Data::Dumper;
use JSON::XS;
use File::Slurp;

use base 'Exporter';

our @EXPORT_OK = qw(reference_fasta_path primer_bedpe_path run_cmds vigor_workflow);

our $ReferenceFasta = "MN908947.fasta";
our $PrimerBedpe = "SC2_200324.bedpe";
our $ArticSchemes = "primer_schemes";
our $VigorWorkflow = "vigor.wf";

sub vigor_workflow
{
    my $mpath = Module::Metadata->find_module_by_name("Bio::P3::SARS2Assembly");
    $mpath =~ s/\.pm$//;

    my $ref = "$mpath/$VigorWorkflow";
    my $txt = read_file($ref);
    my $data = decode_json($txt);
    return $data;
}

sub reference_fasta_path
{
    my $mpath = Module::Metadata->find_module_by_name("Bio::P3::SARS2Assembly");
    $mpath =~ s/\.pm$//;

    my $ref = "$mpath/$ReferenceFasta";
    return $ref;
}

sub primer_bedpe_path
{
    my $mpath = Module::Metadata->find_module_by_name("Bio::P3::SARS2Assembly");
    $mpath =~ s/\.pm$//;

    my $ref = "$mpath/$PrimerBedpe";
    return $ref;
}

sub artic_primer_schemes_path
{
    my $mpath = Module::Metadata->find_module_by_name("Bio::P3::SARS2Assembly");
    $mpath =~ s/\.pm$//;

    my $ref = "$mpath/$ArticSchemes";
    return $ref;
}

sub run_cmds
{
    my(@cmds) = @_;

    print STDERR "Execute:\n";
    for my $c (@cmds)
    {
	if (ref($c) eq 'ARRAY')
	{
	    print STDERR "\t@$c\n";
	}
	elsif(!ref($c))
	{
	    print STDERR "\t$c\n";
	}
    }

    my $ok = run(@cmds);
    $ok or die "Failed running pipeline: \n" . Dumper(\@cmds);
}

1;
