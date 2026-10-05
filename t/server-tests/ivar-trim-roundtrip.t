#
# End-to-end check of the alignment and primer-trimming path used by the onecodex
# recipe, on synthetic reads whose correct output is known exactly.
#
# Mirrors the command chain in service-scripts/sars2-onecodex.pl:
#
#   minimap2 -K 20M -a -x sr            (:262)
#   samtools view -u -h -q 20 -F 4 | samtools sort              (:274)
#   ivar trim -e -q 0 -b <rewritten bed>                        (:291)
#   samtools mpileup --fasta-ref --max-depth 0 --count-orphans
#                    --no-BAQ --min-BQ 0                        (:370)
#   ivar variants -t 0.6                                        (:380)
#
# including the bed rewrite at :135, which keeps only columns 1-4, forces score 60
# and derives strand from the primer NAME rather than the bed's strand column.
#
# Two fixtures:
#   reference  reads cut from the scheme's own reference; nothing should be called
#   variant    reads carrying the five substitutions the v5.4.2 spike-in primers
#              target. Each sits under a primer-binding site, so each must be
#              recovered from the overlapping neighbouring amplicon after trimming.
#
# Also asserts that v5.3.2 and v5.4.2 trim identically. The spike-ins are alternate
# oligos for binding sites that already exist, so they add no new trim interval.
# Compare POS/CIGAR, not whole records: ivar trim writes an XA:i:<n> tag holding the
# index of the matched primer in the bed, which shifts when rows are inserted.
#
# Requires minimap2, samtools and ivar. Skips cleanly if any is missing.
#

use strict;
use warnings;

use Test::More;
use File::Basename;
use File::Spec;
use File::Temp qw(tempdir);
use FindBin;
use lib "$FindBin::Bin/../../lib";

use Bio::P3::SARS2Assembly qw(artic_primer_schemes_path);

my $data     = File::Spec->catdir($FindBin::Bin, '..', 'data');
my $gen      = File::Spec->catfile($data, 'gen_amplicon_reads.py');
my $varfile  = File::Spec->catfile($data, 'v5.4.2-spikein-variants.tsv');
my $root     = artic_primer_schemes_path();
my $scheme   = File::Spec->catdir($root, 'SARS-CoV-2', 'V5.4.2');
my $prev     = File::Spec->catdir($root, 'SARS-CoV-2', 'V5.3.2');

#
# Preconditions.
#
for my $tool (qw(minimap2 samtools ivar))
{
    my $found = `which $tool 2>/dev/null`;
    chomp $found;
    plan skip_all => "$tool is not on PATH; source the runtime's user-env.sh"
        unless $found && -x $found;
}
plan skip_all => "python3 is not on PATH" unless `which python3 2>/dev/null`;
plan skip_all => "V5.4.2 is not in this deployment; run: make artic_schemes"
    unless -d $scheme;

my $bed = File::Spec->catfile($scheme, 'SARS-CoV-2.primer.bed');
my $ref = File::Spec->catfile($scheme, 'SARS-CoV-2.reference.fasta');
plan skip_all => "missing $bed" unless -f $bed;
plan skip_all => "missing $ref" unless -f $ref;

my $dir = tempdir(CLEANUP => 1);

sub run_ok
{
    my($label, @cmd) = @_;
    my $rc = system("@cmd");
    ok($rc == 0, $label) or diag("command failed (rc=$rc): @cmd");
    return $rc == 0;
}

#
# Reproduce the bed rewrite from sars2-onecodex.pl:131-138.
#
sub rewrite_bed
{
    my($in, $out) = @_;
    open(my $i, '<', $in)  or die "$in: $!";
    open(my $o, '>', $out) or die "$out: $!";
    while (<$i>)
    {
        chomp;
        my @x = split(/\t/);
        print $o join("\t", @x[0..3], 60, $x[3] =~ m/LEFT|(F$)/ ? "+" : "-"), "\n";
    }
    close($i); close($o);
    return $out;
}

#
# Align, trim and pile up one read set. Returns the paths it produced.
#
sub pipeline
{
    my($tag, $r1, $r2, $bedfile, $reference) = @_;

    my $trimmed_ref = "$dir/$tag.ref.fa";
    system(qq{perl -pe 's/^(>\\S+).*\$/\\1/' $reference | seqtk trimfq -e 33 - > $trimmed_ref 2>/dev/null})
        == 0 or do {
            # seqtk is optional here; fall back to the untrimmed reference
            system("cp $reference $trimmed_ref");
        };

    run_ok("$tag: minimap2",
           "minimap2 -K 20M -a -x sr -t 4 $trimmed_ref $r1 $r2 -o $dir/$tag.sam 2>/dev/null")
        or return;
    run_ok("$tag: samtools view/sort",
           "samtools view -u -h -q 20 -F 4 $dir/$tag.sam | "
           . "samtools sort --threads 4 -o $dir/$tag.sorted.bam - 2>/dev/null")
        or return;
    system("samtools index $dir/$tag.sorted.bam");

    run_ok("$tag: ivar trim",
           "ivar trim -e -q 0 -i $dir/$tag.sorted.bam -b $bedfile -p $dir/$tag.ivar "
           . "> $dir/$tag.trim.txt 2>&1")
        or return;

    system("samtools sort $dir/$tag.ivar.bam --threads 4 -o $dir/$tag.final.bam 2>/dev/null");
    system("samtools index $dir/$tag.final.bam");
    return "$dir/$tag";
}

#
# ---------------------------------------------------------------- fixture 1: reference
#
my $rbed = rewrite_bed($bed, "$dir/v542.bed");

run_ok("generate reference reads",
       "python3 $gen --ref $ref --bed $bed --r1 $dir/ref_R1.fq --r2 $dir/ref_R2.fq "
       . "--read-len 250 --depth 6 > $dir/gen_ref.log 2>&1");

my $base = pipeline('ref', "$dir/ref_R1.fq", "$dir/ref_R2.fq", $rbed, $ref);

SKIP: {
    skip("reference pipeline did not complete", 4) unless $base && -f "$base.final.bam";

    my $report = do { open my $fh, '<', "$base.trim.txt" or die; local $/; <$fh> };
    like($report, qr/Trimmed primers from 100% /,
         "reference: ivar trimmed a primer from every read")
        or diag($report);
    like($report, qr/0% \(0\) of reads started outside of primer regions/,
         "reference: no read started outside a primer region");

    system("samtools mpileup --fasta-ref $ref --max-depth 0 --count-orphans --no-BAQ "
           . "--min-BQ 0 $base.final.bam > $dir/ref.pileup 2>/dev/null");

    # Every observed base must match the reference. Walk the pileup rather than
    # comparing an ivar consensus positionally: consensus output is not indexed by
    # reference coordinate, so a positional diff against it is meaningless.
    my ($obs, $mismatch) = (0, 0);
    open(my $pu, '<', "$dir/ref.pileup") or die;
    while (<$pu>)
    {
        chomp;
        my @f = split(/\t/);
        next unless @f >= 5;
        my $b = $f[4];
        $b =~ s/\^.//g;            # read start marker plus its mapping quality
        $b =~ s/\$//g;             # read end marker
        $b =~ s/[+-](\d+)/'#' x (length($1) + 1 + $1)/ge;  # blank out indel payloads
        $b =~ s/#//g;
        $obs      += length($b);
        $mismatch += ($b =~ tr/.,//c);
    }
    close($pu);
    cmp_ok($obs, '>', 0, "reference: pileup produced base observations ($obs)");
    is($mismatch, 0, "reference: no pileup base disagrees with the reference");
}

#
# ------------------------------------------------------------------ fixture 2: variant
#
run_ok("generate variant reads",
       "python3 $gen --ref $ref --bed $bed --r1 $dir/var_R1.fq --r2 $dir/var_R2.fq "
       . "--variants $varfile --read-len 250 --depth 6 > $dir/gen_var.log 2>&1");

my $vbase = pipeline('var', "$dir/var_R1.fq", "$dir/var_R2.fq", $rbed, $ref);

SKIP: {
    skip("variant pipeline did not complete", 3) unless $vbase && -f "$vbase.final.bam";

    system("samtools mpileup --fasta-ref $ref --max-depth 0 --count-orphans --no-BAQ "
           . "--min-BQ 0 $vbase.final.bam > $dir/var.pileup 2>/dev/null");
    run_ok("variant: ivar variants",
           "ivar variants -p $dir/var -r $ref -t 0.6 < $dir/var.pileup > /dev/null 2>&1");

    # Expected calls, 1-based, from the 0-based fixture.
    my %want;
    open(my $vf, '<', $varfile) or die "$varfile: $!";
    while (<$vf>)
    {
        next if /^\s*#/ || !/\S/;
        my($pos, $r, $a) = split;
        $want{$pos + 1} = "$r>$a";
    }
    close($vf);

    my %got;
    open(my $tsv, '<', "$dir/var.tsv") or die "no ivar variants output: $!";
    my $hdr = <$tsv>;
    while (<$tsv>)
    {
        my @f = split(/\t/);
        $got{$f[1]} = "$f[2]>$f[3]";
    }
    close($tsv);

    is_deeply(\%got, \%want,
              "variant: exactly the " . scalar(keys %want) . " spike-in target "
              . "substitutions are called, and nothing else")
        or diag("want: " . join(" ", map { "$_=$want{$_}" } sort { $a <=> $b } keys %want)
                . "\ngot:  " . join(" ", map { "$_=$got{$_}" } sort { $a <=> $b } keys %got));

    # The point of the fixture: each of these sits under a primer that ivar clipped,
    # so coverage here can only come from the overlapping neighbouring amplicon.
    my %depth;
    open(my $pu, '<', "$dir/var.pileup") or die;
    while (<$pu>) { my @f = split(/\t/); $depth{$f[1]} = $f[3] if @f >= 4 }
    close($pu);

    my @uncovered = grep { !$depth{$_} } sort { $a <=> $b } keys %want;
    ok(!@uncovered,
       "variant: every under-primer position keeps coverage from a neighbouring amplicon")
        or diag("no depth at: @uncovered");
}

#
# -------------------------------------------------- v5.3.2 / v5.4.2 trimming equivalence
#
SKIP: {
    skip("V5.3.2 not deployed", 1) unless -d $prev;
    my $pbed = File::Spec->catfile($prev, 'SARS-CoV-2.primer.bed');
    skip("V5.3.2 primer bed missing", 1) unless -f $pbed;
    skip("reference pipeline did not complete", 1) unless $base && -f "$base.sorted.bam";

    my $rprev = rewrite_bed($pbed, "$dir/v532.bed");
    system("ivar trim -e -q 0 -i $base.sorted.bam -b $rprev -p $dir/prev.ivar "
           . "> $dir/prev.trim.txt 2>&1");
    system("samtools sort $dir/prev.ivar.bam --threads 4 -o $dir/prev.final.bam 2>/dev/null");

    # POS and CIGAR only: the XA:i: primer-index tag legitimately shifts when the bed
    # gains rows, so comparing whole records would report a spurious difference.
    my $a = `samtools view $dir/prev.final.bam | cut -f1-6 | sort`;
    my $b = `samtools view $base.final.bam     | cut -f1-6 | sort`;
    is($a, $b,
       "v5.3.2 and v5.4.2 trim to identical POS/CIGAR: the spike-ins add no trim interval");
}

done_testing();
