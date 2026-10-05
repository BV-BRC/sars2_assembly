#
# End-to-end confirmation that a primer-site dropout becomes a hole in the
# consensus, and that v5.4.2 fills the holes v5.3.2 leaves.
#
# t/client-tests/primer-dropout-model.t decides which amplicons can still be
# primed on a drifted template; that is the discriminating check and needs no
# tools. This one answers the question that only the real pipeline can: what a
# user actually gets out the other end. The whole onecodex chain is replayed
# (see ivar-trim-roundtrip.t for the per-stage citation), and the assertion is
# on zero-coverage windows and the length of the N run in the consensus.
#
# Measure coverage with `samtools depth -a`, NOT by diffing the consensus against
# the reference positionally: ivar consensus output is not indexed by reference
# coordinate -- it omits the uncovered genome ends -- so an N run reported at
# 3156 in the consensus is the hole at 3234 in the reference. The LENGTH is
# meaningful and is asserted; the offset is not.
#
# The two windows at the genome ends (1-78 and 29841-29903) lie outside the
# outermost amplicons and are present in every run, including the baseline.
#
# Requires minimap2, samtools and ivar. Skips cleanly if any is missing.
#

use strict;
use warnings;

use Test::More;
use File::Spec;
use File::Temp qw(tempdir);
use FindBin;
use lib "$FindBin::Bin/../../lib";

use Bio::P3::SARS2Assembly qw(artic_primer_schemes_path);

my $data      = File::Spec->catdir($FindBin::Bin, '..', 'data');
my $gen       = File::Spec->catfile($data, 'gen_amplicon_reads.py');
my $scenarios = File::Spec->catfile($data, 'dropout-scenarios.tsv');
my $root      = artic_primer_schemes_path();

for my $tool (qw(minimap2 samtools ivar python3))
{
    my $found = `which $tool 2>/dev/null`;
    chomp $found;
    plan skip_all => "$tool is not on PATH; source the runtime's user-env.sh"
        unless $found && -x $found;
}
plan skip_all => "missing $gen"       unless -f $gen;
plan skip_all => "missing $scenarios" unless -f $scenarios;

my %scheme;
for my $v (qw(V5.3.2 V5.4.2))
{
    my $d = File::Spec->catdir($root, 'SARS-CoV-2', $v);
    plan skip_all => "$v is not in this deployment; run: make artic_schemes"
        unless -d $d;
    $scheme{$v} = {
        bed => File::Spec->catfile($d, 'SARS-CoV-2.primer.bed'),
        ref => File::Spec->catfile($d, 'SARS-CoV-2.reference.fasta'),
    };
    plan skip_all => "missing $scheme{$v}{bed}" unless -f $scheme{$v}{bed};
    plan skip_all => "missing $scheme{$v}{ref}" unless -f $scheme{$v}{ref};
}

my $dir = tempdir(CLEANUP => 1);
my $run = 0;

sub sh
{
    my($cmd) = @_;
    my $out = `$cmd 2>&1`;
    die "command failed (rc=$?): $cmd\n$out" if $?;
    return $out;
}

#
# Generate reads under the strict annealing model, push them through the chain,
# and return the zero-coverage windows (1-based inclusive, >=20 bp) plus the
# lengths of the N runs in the consensus.
#
sub pipeline
{
    my($v, $scenario) = @_;
    my $w = File::Spec->catdir($dir, "run" . $run++);
    mkdir($w) or die "cannot mkdir $w: $!";

    my($bed, $ref) = @{$scheme{$v}}{qw(bed ref)};

    sh("python3 '$gen' --ref '$ref' --bed '$bed'"
       . " --scenarios '$scenarios' --scenario '$scenario'"
       . " --anneal-model strict --r1 '$w/r1.fq' --r2 '$w/r2.fq'");

    # The bed rewrite at sars2-onecodex.pl:135 -- cols 1-4, score forced to 60,
    # strand re-derived from the primer NAME rather than taken from the bed.
    open(my $in,  '<', $bed)          or die "cannot read $bed: $!";
    open(my $out, '>', "$w/trim.bed") or die "cannot write $w/trim.bed: $!";
    while (<$in>)
    {
        chomp;
        my @x = split /\t/;
        print $out join("\t", @x[0..3], 60, $x[3] =~ m/LEFT|(F$)/ ? "+" : "-"), "\n";
    }
    close($in);
    close($out);

    sh("minimap2 -K 20M -a -x sr -t 4 '$ref' '$w/r1.fq' '$w/r2.fq' -o '$w/mm.sam'");
    sh("samtools view -u -h -q 20 -F 4 '$w/mm.sam'"
       . " | samtools sort --threads 4 -o '$w/s.bam' -");
    sh("samtools index '$w/s.bam'");
    sh("ivar trim -e -q 0 -i '$w/s.bam' -b '$w/trim.bed' -p '$w/t'");
    sh("samtools sort --threads 4 -o '$w/is.bam' '$w/t.bam'");
    sh("samtools index '$w/is.bam'");
    sh("samtools mpileup --fasta-ref '$ref' --max-depth 0 --count-orphans"
       . " --no-BAQ --min-BQ 0 '$w/is.bam' > '$w/pileup'");
    sh("ivar consensus -p '$w/c' -m 3 -t 0.6 -n N < '$w/pileup'");

    # Zero-coverage windows, collapsed into runs. depth -a emits every reference
    # position, so these ARE reference coordinates.
    my @zero;
    open(my $d, '-|', "samtools depth -a '$w/is.bam'") or die "samtools depth: $!";
    while (<$d>)
    {
        my(undef, $pos, $depth) = split;
        next if $depth != 0;
        if (@zero && $pos == $zero[-1][1] + 1) { $zero[-1][1] = $pos }
        else                                   { push(@zero, [$pos, $pos]) }
    }
    close($d);
    my @windows = map { "$_->[0]-$_->[1]" }
                  grep { $_->[1] - $_->[0] + 1 >= 20 } @zero;

    open(my $c, '<', "$w/c.fa") or die "cannot read $w/c.fa: $!";
    my $cons = join('', map { chomp; /^>/ ? () : $_ } <$c>);
    close($c);
    my @nruns = map { length } grep { length($_) >= 20 } ($cons =~ /(N+)/g);

    return { windows => \@windows, nruns => \@nruns };
}

#
# The ends of the genome fall outside the outermost amplicons and are never
# covered by any scheme. Everything below is stated relative to this.
#
my @ENDS = ('1-78', '29841-29903');

#
# Baseline: undrifted reference. Both schemes cover the genome identically, so
# adding five oligos costs nothing when there is no drift to answer.
#
for my $v (qw(V5.3.2 V5.4.2))
{
    my $r = pipeline($v, 'none');
    is_deeply($r->{windows}, \@ENDS, "$v/none: no coverage hole beyond the genome ends");
    is_deeply($r->{nruns},   [],     "$v/none: consensus has no N run");
}

#
# One drift scenario per affected amplicon. v5.3.2 loses the amplicon and the
# hole is exactly the span its neighbours' inserts do not reach; v5.4.2 has a
# spiked-in oligo that still primes, and the hole closes completely.
#
# The expected windows are the tiling geometry: dropping amplicon N leaves
# prev.insert_end .. next.insert_start uncovered. They are spelled out rather
# than computed so that a scheme edit which changed the tiling would have to be
# acknowledged here.
#
my @CASES = (
    { scenario => 'amp11',  amp => 11, hole => '3234-3540',  len => 307 },
    { scenario => 'amp69a', amp => 69, hole => '21373-21607', len => 235 },
    { scenario => 'amp70',  amp => 70, hole => '21697-21894', len => 198 },
    { scenario => 'amp74',  amp => 74, hole => '22840-23109', len => 270 },
);

for my $c (@CASES)
{
    my($s, $amp, $hole, $len) = @{$c}{qw(scenario amp hole len)};

    my $old = pipeline('V5.3.2', $s);
    is_deeply($old->{windows}, [$ENDS[0], $hole, $ENDS[1]],
              "$s: v5.3.2 loses amplicon $amp -- $len bp of zero coverage at $hole");
    is_deeply($old->{nruns}, [$len],
              "$s: v5.3.2 consensus has one $len bp N run");

    my $new = pipeline('V5.4.2', $s);
    is_deeply($new->{windows}, \@ENDS,
              "$s: v5.4.2 keeps amplicon $amp -- the $hole hole is gone");
    is_deeply($new->{nruns}, [],
              "$s: v5.4.2 consensus has no N run");
}

done_testing();
