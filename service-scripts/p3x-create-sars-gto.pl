#
# Create a GTO from a sars2 assembly directory
#

use strict;
use Data::Dumper;
use Getopt::Long::Descriptive;
use File::Slurp;
use JSON::XS;
use GenomeTypeObject;
use gjoseqlib;

my($opt, $usage) = describe_options("%c %o contigs sra-metadata-json output.gto",
				    ["help" => "Show this help message."],
				   );
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 3;

my $contigs_file = shift;
my $metadata_file = shift;
my $output = shift;

my @contigs;
open(CTG, "<", $contigs_file) or die "Cannot open $contigs_file: $!";
while (my($id, $def, $seq) = read_next_fasta_seq(\*CTG))
{
    push(@contigs, {
	id => $id,
	dna => $seq,
    });
}

my $metadata = eval { decode_json(scalar read_file($metadata_file)); };
$metadata or die "Could not read and parse metadata file $metadata_file: $@";

ref($metadata) eq 'ARRAY' or die "Unexpected non-array metadata\n";
$metadata = $metadata->[0];

my $gto = GenomeTypeObject->new;
$gto->add_contigs(\@contigs);

$gto->{id} = $metadata->{accession};
$gto->{sra_metadata} = $metadata;
$gto->{ncbi_taxonomy_id} = $metadata->{sample_taxon} // 2;
$gto->{scientific_name} = $metadata->{sample_organism} // "Unknown sp.";
$gto->{domain} = 'Viruses';
$gto->{genetic_code} = 1;

$gto->destroy_to_file($output);
    

    
