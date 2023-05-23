#
# Reformat primers to our convention.
#
# Usage: rewrite-primers.pl scheme-dir
#
#
# This replaces the following confusing and fragile make code:
#	cd lib/Bio/P3/SARS2Assembly; \
#	for s in primer_schemes/nCoV-2019/V*; do \
#	   (cd $$s; \
#	   perl -ne 'my @x=split m/\t/; print join("\t",@x[0..3], 60, $x[3]=~m/LEFT/?"+":"-"),"\n";' \
#		?(nCoV-2019|SARS-CoV-2).scheme.bed) > ARTIC-`basename $$s`.bed; \
#	done
#

use strict;
use File::Basename;

@ARGV == 2 or die "Usage: $0 scheme-dir dest-dir";

my $scheme_dir = shift;
my $dest_dir = shift;

-d $scheme_dir or die "Scheme dir $scheme_dir does not exist";
-d $dest_dir or die "Destination dir $dest_dir does not exist";

for my $vers_dir (<$scheme_dir/V*>)
{
    my $vers = basename($vers_dir);
    my $bed = "$vers_dir/nCoV-2019.scheme.bed";
    if (! -f $bed)
    {
	$bed = "$vers_dir/SARS-CoV-2.scheme.bed";
    }
    if (! -f $bed)
    {
	die "Could not find bed file in $vers_dir\n";
    }
    my $bed_out = "$dest_dir/ARTIC-$vers.bed";
    open(BED_IN, "<", $bed) or die "Cannot open $bed: $!";
    open(BED_OUT, ">", $bed_out) or die "Cannot open $bed_out: $!";
    print "Rewrite $bed to $bed_out\n";
    while (<BED_IN>)
    {
	my @x = split m/\t/;
	print BED_OUT join("\t", @x[0..3], 60, $x[3] =~ m/LEFT/ ? "+" : "-"), "\n";
    }
    close(BED_IN);
    close(BED_OUT);
}
