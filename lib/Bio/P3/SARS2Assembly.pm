package Bio::P3::SARS2Assembly;

#
# Mostly empty package; used as a placeholder 
#

use strict;
use Module::Metadata;
use IPC::Run qw(run timeout);
use Data::Dumper;
use JSON::XS;
use Time::HiRes 'gettimeofday';
use File::Slurp;

use base 'Exporter';

our @EXPORT_OK = qw(reference_fasta_path primer_bedpe_path run_cmds vigor_workflow report_template
		    reference_gff_path reference_spike_aa_path mpath
		    add_variants_to_gto add_quality_estimate_to_gto
		    artic_bed artic_reference);

our $ReferenceSpikeAA = "YP_009724390.1.aa.fa";
our $ReferenceFasta = "MN908947.fasta";
our $ReferenceGFF = "GCF_009858895.2_ASM985889v3_genomic.gff";
our $PrimerBedpe = "SC2_200324.bedpe";
our $ArticSchemes = "primer_schemes";
our $VigorWorkflow = "vigor.wf";
our $ReportTemplate = "report.tt";

sub mpath
{
    my $mpath = Module::Metadata->find_module_by_name("Bio::P3::SARS2Assembly");
    $mpath =~ s/\.pm$//;
    return $mpath;
}

sub vigor_workflow
{
    my $ref = mpath() . "/$VigorWorkflow";
    my $txt = read_file($ref);
    my $data = decode_json($txt);
    return $data;
}

sub reference_fasta_path
{
    my $ref = mpath() . "/$ReferenceFasta";
    return $ref;
}

sub reference_spike_aa_path
{
    my $ref = mpath() . "/$ReferenceSpikeAA";
    return $ref;
}

sub reference_gff_path
{
    my $ref = mpath() . "/$ReferenceGFF";
    return $ref;
}

sub report_template
{
    my $ref = mpath() . "/$ReportTemplate";
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

sub artic_reference
{
    my($vers) = @_;
    my $ref = mpath() . "/$ArticSchemes/nCoV-2019/V$vers/nCoV-2019.reference.fasta";
    return $ref;
}

sub artic_bed
{
    my($vers) = @_;
    my $ref = mpath() . "/ARTIC-V$vers.bed";
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

    my $start = gettimeofday;
    my $ok = run(@cmds);
    my $end = gettimeofday;
    my $elap = $end - $start;
    print STDERR "Run returns $ok $? elapsed=$elap\n";
    $ok or die "Failed running pipeline: \n" . Dumper(\@cmds);
}

sub add_variants_to_gto
{
    my($file, $gto) = @_;
    #
    # Read reference.
    #

    my %gene;

    my @prots;
    open(R, "<", reference_gff_path) or die "Cannot open reference " . reference_gff_path . ": $!";
    while (<R>)
    {
	next if /^#/;
	    chomp;
	my($contig, $src, $type, $start, $end, $score, $strand, $frame, $attr) = split(/\t/);
	next unless $type eq 'CDS';

	my $attrh = {};
	for my $att (split(/;/, $attr))
	{
	    my($key, $val) = split(/=/, $att, 2);
	    $attrh->{$key} = $val;
	}

	$gene{$attrh->{ID}} = $attrh->{gene};
	push(@prots, [$start, $end, $attrh]);
    }
    close(R);
    
    open(V, "<", $file) or die "Cannot open variants file $file: $!\n";
    my $hdr = <V>;
    chomp $hdr;
    my @keys = split(/\t/, $hdr);
    my %idx;
    $idx{$keys[$_]} = $_ foreach (0..$#keys);

    my %snplist;
    while (<V>)
    {
	chomp;
	my @dat = split(/\t/);
	my $pos = 0 + $dat[$idx{POS}];
	my $freq = 0 + $dat[$idx{ALT_FREQ}];
	my $feat = $dat[$idx{GFF_FEATURE}];
	my $ref_aa = $dat[$idx{REF_AA}];
	my $alt_aa = $dat[$idx{ALT_AA}];
	my $ref = $dat[$idx{REF}];
	my $alt = $dat[$idx{ALT}];

	my $snp = {
	    pos => $pos,
	    ref => $ref,
	    alt => $alt,
	    freq => $freq,
	};

	my $prot;
	for my $item (@prots)
	{
	    my($s, $e, $attr) = @$item;
	    if ($pos >= $s && $pos < $e)
	    {
		$prot = $item;
		last;
	    }
	}

	if ($prot)
	{
	    my($s, $e, $attr) = @$prot;

	    #
	    # If we have an AA reference that is nonsynonomous, record it
	    #
	    if ($ref_aa)
	    {
		$snp->{ref_aa} = $ref_aa;
		$snp->{alt_aa} = $alt_aa;
		
		my $off1 = $pos - $s;
		my $off = int($off1 / 3) + 1;

		$snp->{feature_pos} = $off;
	    }

	    push(@{$snplist{$attr->{ID}}}, $snp);
	}
	else
	{
	    #
	    # Record the dna sub as intergenic.
	    #
	    push(@{$snplist{Intergenic}}, $snp);
	}
    }
    close(V);
    my $vlist = [];
    for my $id (sort keys %snplist)
    {
	my $gene = $gene{$id};
	my $tv = {
	    reference => $id,
	    gene => $gene,
	    snps => $snplist{$id} // [],
	};
	push(@$vlist, $tv);
    }
    push(@{$gto->{computed_variants}}, {
	tool => "assembly pipeline",
	variants => $vlist,
    });
}

sub add_quality_estimate_to_gto
{
    my($gto) = @_;

    my $max_n = 0;
    my $total_nc = 0;
    my $total_contigs = 0;
    for my $ctg ($gto->contigs)
    {
	$total_contigs += length($ctg->{dna});
	while ($ctg->{dna} =~ /([Nn]+)/g)
	{
	    my $nc = length($1);
	    $total_nc += $nc;
	    $max_n = $nc if $nc > $max_n;
	}
    }

    $gto->{quality}->{contig_ambig_count} = $total_nc;
    $gto->{quality}->{contig_ambig_fraction} = $total_nc / $total_contigs;
    $gto->{quality}->{contig_longest_ambig_run} = $max_n;
}

1;


