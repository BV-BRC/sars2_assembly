
#
# Given a file with a list of SRA IDs and a set of base directories for output,
# search for status.
#

use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;

my($opt, $usage) = describe_options("%c %o sra-ids base-dir [base-dir...]",
				    ["help|h" => "Show this help message."],
				   );
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV < 2;

my $sra_id_file = shift;
my @base_dirs = @ARGV;

my %ids;
my @ids;
my %prefixes;
open(S, "<", $sra_id_file) or die "Cannot read id file $sra_id_file: $!";

while (<S>)
{
    chomp;
    my($x) = split(/\t/);
    if ($x =~ /^([A-Z]{3}\d{4})\d{3}/)
    {
	my $id = $x;
	my $prefix = $1;
	$prefixes{$prefix} = 1;
	push(@ids, $id);
	$ids{$id} = 1;
    }
}
close(S);

my %to_process = %ids;
for my $base (@base_dirs)
{
    if (!opendir(B, $base))
    {
	warn "Failure opening directory $base: $!";
	next;
    }
    while (my $p = readdir(B))
    {
	next unless $prefixes{$p};
	my $prefix_dir = "$base/$p";
	if (!opendir(P, $prefix_dir))
	{
	    warn "Failure opening dir $prefix_dir: $!";
	    next;
	}
	while (my $sra = readdir(P))
	{
	    next unless $ids{$sra};
	    if (-s "$prefix_dir/$sra/$sra.gto")
	    {
		delete $to_process{$sra};
	    }
	}
    }
}

for my $id ( grep { $to_process{$_} } @ids)
{
    print "$id\n";
}
