#
# Does v5.4.2 actually fix the amplicon dropouts it was released for?
#
# ARTIC v5.4.2 spikes five alternate oligos into the v5.3.2 400 bp pool to restore
# amplification against JN.1-era drift. The five substitutions they answer to sit
# UNDER primer-binding sites, so the failure mode is not a miscalled base -- it is
# the whole amplicon failing to amplify.
#
# This is the discriminating half of the dropout check and needs no bioinformatics
# tools: it asks, for each drifted template, whether each scheme still has an oligo
# that can prime. t/server-tests/primer-dropout.t then confirms that the predicted
# dropouts really do become holes in the consensus.
#
# The annealing rule is a model, and it is the one assumption here that is not a
# measured fact. Two are exercised:
#
#   strict       any primer/template mismatch kills the oligo. Upper bound on
#                sensitivity; treats all five spike-ins as answering a real dropout.
#   three-prime  a mismatch within 5 nt of the 3' end blocks extension, and >1
#                mismatch anywhere kills the oligo; a lone mismatch further 5' is
#                tolerated. Closer to the chemistry, and deliberately stricter about
#                what counts as a dropout -- under it only amplicons 11 and 74
#                qualify, because only their substitutions are 3'-proximal. The
#                other three depress Tm rather than blocking extension, which a
#                binary in-silico model cannot represent faithfully.
#
# Both agree on the result that matters: v5.4.2's dropout set is a subset of
# v5.3.2's in every scenario, and empty in every single-mutation one.
#
# An amplicon amplifies when ANY of its LEFT and ANY of its RIGHT oligos anneal --
# that is the whole mechanism by which a spiked-in alternate rescues a dropout, and
# why these are additions to the pool rather than replacements.
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

plan skip_all => "python3 is not on PATH" unless `which python3 2>/dev/null`;
plan skip_all => "missing $gen"           unless -f $gen;
plan skip_all => "missing $scenarios"     unless -f $scenarios;

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

#
# Run the generator for one scheme/scenario/model and return the sorted list of
# amplicons that failed to amplify. Reads go to /dev/null -- only the annealing
# decision is under test here.
#
sub dropouts
{
    my($v, $scenario, $model) = @_;
    my $report = File::Spec->catfile($dir, "report.tsv");
    my @cmd = ('python3', $gen,
               '--ref'            => $scheme{$v}{ref},
               '--bed'            => $scheme{$v}{bed},
               '--scenarios'      => $scenarios,
               '--scenario'       => $scenario,
               '--anneal-model'   => $model,
               '--dropout-report' => $report,
               '--r1'             => '/dev/null',
               '--r2'             => '/dev/null');
    my $out = `@cmd 2>&1`;
    die "generator failed for $v/$scenario/$model (rc=$?):\n$out" if $?;

    open(my $fh, '<', $report) or die "cannot read $report: $!";
    my @drop;
    while (<$fh>)
    {
        next if $. == 1;
        my($amp, $status) = split /\t/;
        push(@drop, $amp) if $status eq 'DROPOUT';
    }
    close($fh);
    return [sort { $a <=> $b } @drop];
}

#
# Baseline: on the undrifted reference neither scheme loses anything, under either
# model. This is what says the five extra oligos introduce no regression of their
# own -- a spike-in that mismatched the reference everywhere would show up here.
#
for my $model (qw(strict three-prime))
{
    for my $v (qw(V5.3.2 V5.4.2))
    {
        is_deeply(dropouts($v, 'none', $model), [],
                  "$v/$model: no dropouts on the undrifted reference");
    }
}

#
# The five drift scenarios, one per independent substitution. Under the strict
# model each one knocks out exactly the amplicon whose primer it sits under in
# v5.3.2, and v5.4.2 rescues all five.
#
my %expect_strict = (
    amp11  => 11,
    amp69a => 69,
    amp69b => 69,
    amp70  => 70,
    amp74  => 74,
);

for my $scenario (sort keys %expect_strict)
{
    my $amp = $expect_strict{$scenario};
    is_deeply(dropouts('V5.3.2', $scenario, 'strict'), [$amp],
              "strict/$scenario: v5.3.2 loses amplicon $amp");
    is_deeply(dropouts('V5.4.2', $scenario, 'strict'), [],
              "strict/$scenario: v5.4.2 rescues amplicon $amp");
}

#
# Under the mechanistic model only the two 3'-proximal substitutions qualify as
# dropouts. Pinning the negative cases too: if a future scheme edit moved one of
# the other three closer to a 3' end, that is a change in kind and should be
# noticed rather than silently absorbed.
#
my %expect_three_prime = (
    amp11  => [11],
    amp69a => [],
    amp69b => [],
    amp70  => [],
    amp74  => [74],
);

for my $scenario (sort keys %expect_three_prime)
{
    my $want = $expect_three_prime{$scenario};
    is_deeply(dropouts('V5.3.2', $scenario, 'three-prime'), $want,
              "three-prime/$scenario: v5.3.2 drops " .
              (@$want ? "amplicon $want->[0]" : "nothing (mismatch is not 3'-proximal)"));
    is_deeply(dropouts('V5.4.2', $scenario, 'three-prime'), [],
              "three-prime/$scenario: v5.4.2 drops nothing");
}

#
# v5.4.2 never does worse than v5.3.2. This is the claim that does not depend on
# which annealing model you believe, and it is the one that justifies promoting
# v5.4.2 to the default.
#
for my $model (qw(strict three-prime))
{
    for my $scenario (qw(none amp11 amp69a amp69b amp70 amp74 all))
    {
        my %old = map { $_ => 1 } @{ dropouts('V5.3.2', $scenario, $model) };
        my @worse = grep { !$old{$_} } @{ dropouts('V5.4.2', $scenario, $model) };
        is_deeply(\@worse, [],
                  "$model/$scenario: v5.4.2 drops nothing that v5.3.2 kept");
    }
}

#
# The all-at-once genome is a chimera and is pinned as such.
#
# Amplicon 69's two spike-ins are alternates for two DIFFERENT single mutations:
# _69_RIGHT_1 answers 21710, _69_RIGHT_2 answers 21717, and neither carries both.
# So a template with both substitutions matches no oligo in either scheme and
# amplicon 69 is legitimately lost -- an artifact of a genome no real lineage
# carries, not a defect in v5.4.2. Asserted rather than avoided: if ARTIC ever
# ships a double-substitution oligo, this test is where that shows up.
#
is_deeply(dropouts('V5.3.2', 'all', 'strict'), [11, 69, 70, 74],
          "strict/all: v5.3.2 loses all four drifted amplicons");
is_deeply(dropouts('V5.4.2', 'all', 'strict'), [69],
          "strict/all: v5.4.2 rescues 11/70/74; 69 needs an oligo carrying both " .
          "21710 and 21717, which no lineage requires and ARTIC did not ship");

done_testing();
