#
# Primer scheme integration checks.
#
# These cover the seam between the primer-schemes data and the code that consumes
# it, using the real path-resolution helpers so the lookup itself is exercised.
# Deep validation of the scheme data (md5s, primer sequence vs reference
# coordinate, amplicon topology) lives at the source, in the primer-schemes
# repository: tests/validate_schemes.py.
#
# This runs against the DEPLOYED copy under lib/Bio/P3/SARS2Assembly/primer_schemes,
# which is what sars2-onecodex.pl reads. If that tree is absent or predates a
# scheme, the affected tests skip with a pointer to `make artic_schemes` rather
# than reporting a data error.
#

use strict;
use warnings;

use Test::More;
use JSON::XS;
use File::Basename;
use File::Spec;
use FindBin;
use lib "$FindBin::Bin/../../lib";

use Bio::P3::SARS2Assembly qw(manifest artic_primer_schemes_path);

my $manifest_path = manifest();
my $schemes_root  = artic_primer_schemes_path();

unless (-f $manifest_path)
{
    plan skip_all => "primer schemes are not deployed (no $manifest_path); run: make artic_schemes";
}

my $raw = do { open my $fh, '<', $manifest_path or die "$manifest_path: $!"; local $/; <$fh> };

my $manifest = eval { decode_json($raw) };
ok(defined $manifest, "bvbrc_manifest.json is valid JSON")
    or BAIL_OUT("cannot parse $manifest_path: $@");

my $org = $manifest->{organisms}[0];
ok(defined $org && $org->{name} eq 'SARS-CoV-2', "manifest declares the SARS-CoV-2 organism");

my $primers = $org->{primers};
ok(ref $primers eq 'HASH' && keys %$primers, "manifest declares primer kits");

#
# Duplicate keys parse as last-wins, silently masking earlier bodies. The manifest
# carried five triplicated kits for years; one of the masked bodies pointed at
# filenames from a different kit. Guard against a reintroduction.
#
{
    my @dup;
    for my $kit (sort keys %$primers)
    {
        my $n = () = $raw =~ /"\Q$kit\E"\s*:/g;
        push(@dup, "$kit (x$n)") if $n > 1;
    }
    ok(!@dup, "no primer kit key is declared more than once")
        or diag("duplicated keys: @dup\n"
                . "this is fixed upstream in BV-BRC/primer-schemes; if the deployed\n"
                . "tree predates that fix, refresh it:\n"
                . "  rm -rf lib/Bio/P3/SARS2Assembly/primer_schemes && make artic_schemes");
}

#
# Every file the manifest promises must actually be there. sars2-onecodex.pl dies
# at run time otherwise, which is a bad way to discover a typo'd or stale entry.
#
{
    my @missing;
    for my $kit (sort keys %$primers)
    {
        my $body = $primers->{$kit};
        for my $scheme (@{$body->{schemes}})
        {
            for my $field (qw(primers reference))
            {
                my $path = File::Spec->catfile($schemes_root, $body->{path},
                                               $scheme->{version}, $scheme->{$field});
                push(@missing, "$kit/$scheme->{version}: $path") unless -f $path;
            }
        }
    }
    ok(!@missing, "every path declared in the manifest exists on disk")
        or diag("missing:\n  " . join("\n  ", @missing));
}

#
# sars2-onecodex.pl:135 does not read the strand column. It derives strand from the
# primer NAME:
#
#     $x[3] =~ m/LEFT|(F$)/ ? "+" : "-"
#
# so a scheme whose names defeat that regex would have its reads trimmed on the
# wrong end, silently. Two naming families are in play: ARTIC-style _LEFT/_RIGHT,
# and swift-style ...F/...R. Check the derivation against the bed's own strand
# column wherever the bed has one.
#
{
    my $checked = 0;
    for my $kit (sort keys %$primers)
    {
        my $body = $primers->{$kit};
        for my $scheme (@{$body->{schemes}})
        {
            my $bed = File::Spec->catfile($schemes_root, $body->{path},
                                          $scheme->{version}, $scheme->{primers});
            next unless -f $bed;

            open(my $fh, '<', $bed) or die "$bed: $!";
            my (@wrong, $rows_with_strand);
            while (<$fh>)
            {
                s/\r?\n\z//;
                next unless length;
                next if /^#/;
                my @f = split(/\t/, $_, -1);
                next unless @f >= 6 && defined $f[5] && $f[5] =~ /^[+-]$/;
                $rows_with_strand++;
                my $derived = $f[3] =~ m/LEFT|(F$)/ ? '+' : '-';
                push(@wrong, "$f[3] bed=$f[5] derived=$derived") if $derived ne $f[5];
            }
            close($fh);

            next unless $rows_with_strand;
            $checked++;
            ok(!@wrong,
               "$kit/$scheme->{version}: /LEFT|(F\$)/ reproduces the bed strand column"
               . " ($rows_with_strand primers)")
                or diag(join("\n  ", @wrong[0 .. ($#wrong > 4 ? 4 : $#wrong)]));
        }
    }
    cmp_ok($checked, '>', 0, "at least one scheme had a strand column to check against");
}

#
# swift has no strand column, and is the reason the regex carries the (F$) branch
# at all. Assert the derivation at least partitions it the way the names imply.
#
SKIP: {
    my $body = $primers->{swift} or skip("no swift kit in manifest", 1);
    my $scheme = $body->{schemes}[0] or skip("no swift scheme", 1);
    my $bed = File::Spec->catfile($schemes_root, $body->{path},
                                  $scheme->{version}, $scheme->{primers});
    skip("swift bed not deployed", 1) unless -f $bed;

    my ($fwd, $rev, $other) = (0, 0, 0);
    open(my $fh, '<', $bed) or die "$bed: $!";
    while (<$fh>)
    {
        s/\r?\n\z//;
        next unless length;
        my @f = split(/\t/, $_, -1);
        my $derived = $f[3] =~ m/LEFT|(F$)/ ? '+' : '-';
        if    ($f[3] =~ /F$/) { $derived eq '+' ? $fwd++ : $other++ }
        elsif ($f[3] =~ /R$/) { $derived eq '-' ? $rev++ : $other++ }
        else                  { $other++ }
    }
    close($fh);
    is($other, 0, "swift: every primer name ends in F or R and derives the matching strand")
        or diag("misclassified or unrecognised: $other");
    cmp_ok($fwd, '==', $rev, "swift: forward and reverse primer counts match ($fwd/$rev)");
}

#
# ARTIC version registration. sars2-onecodex.pl:113 takes $schemes->[-1] when the
# caller passes no primer_version, so order is what selects the default.
#
{
    my $artic = $primers->{ARTIC};
    ok(defined $artic, "ARTIC kit is registered");

    my @versions = map { $_->{version} } @{$artic->{schemes}};
    cmp_ok(scalar @versions, '>', 0, "ARTIC declares at least one version (@versions)");

    my %seen;
    my @repeated = grep { $seen{$_}++ } @versions;
    ok(!@repeated, "no ARTIC version is declared twice")
        or diag("repeated: @repeated");

  SKIP: {
        my $deployed = File::Spec->catdir($schemes_root, $artic->{path}, 'V5.4.2');
        skip("V5.4.2 is not in this deployment; run: make artic_schemes", 2)
            unless -d $deployed;

        ok(grep({ $_ eq 'V5.4.2' } @versions), "V5.4.2 is registered in the manifest");
        is($versions[-1], 'V5.4.2',
           "V5.4.2 is last, so it is the default for ARTIC jobs with no primer_version")
            or diag("order: @versions");
    }
}

done_testing();
