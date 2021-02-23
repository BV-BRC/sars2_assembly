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
use Bio::P3::SARS2Assembly qw(reference_gff_path add_variants_to_gto);

my($opt, $usage) = describe_options("%c %o contigs sra-metadata-json output.gto",
				    ["variants=s" => "Variants vcf file"],
				    ["accession=s" => "SRA accession number"],
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
if (!$metadata)
{
    warn "Could not read and parse metadata file $metadata_file: $@";
    $metadata = {
	(defined($opt->accession) ? (run_id => $opt->accession, accession => $opt->accession) : ()),
	sample_taxon => 2697049,
	sample_organism => "Severe acute respiratory syndrome coronavirus 2",
    };
    $metadata = [$metadata];
}

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

#
# Process variants if provided
#

if ($opt->variants)
{
    add_variants_to_gto($opt->variants, $gto);
}


$gto->destroy_to_file($output, { canonical => 1 });
    

    
