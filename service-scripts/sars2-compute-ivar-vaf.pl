
#
# Given an ivar variants file on stdin, compute VAF as
# (ALT_DP + ALT_RP) / (ALT_DP + ALT_RP + REF_DP + REF_RV)
# and add a column.
#

use strict;
use Data::Dumper;

my $hdrs = <>;
chomp $hdrs;
my @hdrs = split(/\t/, $hdrs);
my %idx;
for my $i (0..$#hdrs)
{
    $idx{$hdrs[$i]} = $i;
}

my $i_alt_dp = $idx{ALT_DP};
my $i_alt_rp = $idx{ALT_RP};

my $i_ref_dp = $idx{REF_DP};
my $i_ref_rp = $idx{REF_RP};

print join("\t", @hdrs, 'VAF'), "\n";

while (<>)
{
    chomp;
    my @dat = split(/\t/);
    my $s = $dat[$i_alt_dp] + $dat[$i_alt_rp];
    my $vaf = $s / ($s + $dat[$i_ref_dp] + $dat[$i_ref_rp]);
    print join("\t", @dat, sprintf("%.4f", $vaf)), "\n";
}
