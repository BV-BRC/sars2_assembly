
=head1 NAME

    sars2-onecodex - OneCodex protocol for SARS2 assembly 
    
=head1 SYNOPSIS

    sars2-onecodex PE-read1 PE-read2 output-base output-dir
    sars2-onecodex SE-read output-base output-dir

=head1 DESCRIPTION

Assemble the input reads using the OneCodex protocol
L<https://github.com/onecodex/sars-cov-2>.    

Adapted from https://github.com/onecodex/sars-cov-2/blob/master/covid19_call_variants.sh

The generated FASTA consensus sequence is written to output-dir/output-base.fasta.

=cut

use strict;
use Getopt::Long::Descriptive;
use IPC::Run qw(run timeout start);
use Bio::P3::SARS2Assembly qw(manifest artic_reference artic_bed run_cmds reference_gff_path artic_primer_schemes_path);
use JSON::XS;
use Data::Dumper;
use File::Basename;
use File::Temp;
use File::Slurp;
use Time::HiRes 'gettimeofday';
use gjoseqlib;
use Bio::P3::CmdRunner;
use PDL;
use PDL::Stats::Basic;
use PDL::Ufunc;

$ENV{PATH} = "$ENV{KB_RUNTIME}/samtools-1.11/bin:$ENV{KB_RUNTIME}/bcftools-1.9/bin:$ENV{PATH}";

#
# The manifest file holds the definitions of the available organisms and primer sets. Use it to
# validate the primer chosen by the user (and to list the available primer sets).
#

my $manifest = decode_json(scalar read_file(manifest));
my($primer_info) = grep { $_->{name} = "SARS-CoV-2" } @{$manifest->{organisms}};
my $primers = $primer_info->{primers};
my @primer_names = sort keys %$primers;

my($opt, $usage) = describe_options("%c %o output-base output-dir",
				    ['pe-read-1|1=s@' => "Paired-end mate 1 file", { default => [] }],
				    ['pe-read-2|2=s@' => "Paired-end mate 2 file", { default => [] }],
				    ['se-read|U=s@' => "Single-end read file", { default => [] }],
				    ["output-name|n=s" => "Output name for sequence (in the fasta file). Defaults to output-base"],
				    ["threads|j=i" => "Number of threads to use", { default => 1 }],
				    ["min-quality|q=i" => "Minimum read quality", { default => 20 }],
				    ["max-depth|d=i" => "Maxmium depth to use in mpileup", { default => 0 }],
				    ["min-depth|D=i" => "Minimum depth for consensus base call", { default => 3 }],
				    ["bed-file=s" => "Use the given primer BED file. The --reference parameter must also be specified"],
				    ["reference=s" => "Use the given reference FASTA file. THe --bed-file parameter must also be specified"],
				    ["primers=s" => "Use these primers. Choices are @primer_names"],
				    ["primer-version=s" => "Use the specfiied version of the chosen primers. Default is the latest available version"],
				    ["length-threshold|l=i" => "Max length to be considered short read sequencing", { default => 600 }],
				    ["nanopore" => "Force use of nanopore mapping method"],
				    ["keep-intermediates|k" => "Save all intermediate files"],
				    ["delete-reads" => "Delete reads when they have been processed"],
				    ["samtools-sort-timeout=i" => "Timeout for samtools sort", { default => 960 }],
				    ["help|h"      => "Show this help message"],
				    );


print($usage->text), exit 0 if $opt->help;
die($usage->text) unless @ARGV == 2;

$opt->primers or die "Primers must be defined using the --primers flag. Available primers are @primer_names\n";

my($reference, $bed_file);
if ($opt->bed_file || $opt->reference)
{
    if (!$opt->bed_file || !$opt->reference)
    {
	die "If either of --bed-file or --reference is specfied, both must be specified\n";
    }
    if ($opt->primers)
    {
	die "Primers may be specified by only one of --bed-file and --primers\n";
    }

    $reference = $opt->reference;
    $bed_file = $opt->bed_file;
}
else
{
    #
    # Find our primer data.
    #
    my $primer = $primers->{$opt->primers};
    $primer or die "Chosen primers " . $opt->primers . " not available\n";
    
    my $scheme;
    my $schemes = $primer->{schemes};
    if ($opt->primer_version)
    {
	($scheme) = grep { $_->{version} eq $opt->primer_version } @$schemes;
	if (!$scheme)
	{
	    my @avail = map { $_->{version} } @$schemes;
	    die "Version " . $opt->primer_version . " not available for primers " . $opt->primers . ". Available versions: @avail\n";
	}
    }
    else
    {
	$scheme = $schemes->[-1];
    }
    
    my $path = artic_primer_schemes_path . "/$primer->{path}/$scheme->{version}";
    $reference = "$path/$scheme->{reference}";
    $bed_file = "$path/$scheme->{primers}";

}

-f $reference or die "Cannot read reference $reference\n";

print STDERR "Processing bed file $bed_file\n";
if (! -f $bed_file)
{
    die "Bed file $bed_file is missing\n";
}
my $bed_tmp = File::Temp->new;
open(I, "<", $bed_file) or die "cannot read $bed_file: $!\n";
while (<I>)
{
    my @x = split m/\t/;
    my $l = join("\t", @x[0..3], 60, $x[3]=~m/LEFT/ ? "+" : "-") . "\n";
    print $bed_tmp $l;
    print STDERR $l;
}
close(I);
close($bed_tmp);

my @stats;

push(@stats,
     [reference => $reference],
     [primers => $bed_file],
     );

my $runner = Bio::P3::CmdRunner->new;

my @se_read_files = @{$opt->se_read};

#
# Unpack the lists of -1 and -2 reads into pairs and verify the counts match.
#
# We also construct the @inputs list which is just the linearized list of
# files. We assume the pair files are named such that minimap2 is happy with them.
#

my @inputs = @se_read_files;

my @pe_read_files;

my $n1 = @{$opt->pe_read_1};
my $n2 = @{$opt->pe_read_2};

if ($n1 != $n2)
{
    die "Mismatch in paired end read file counts ($n1 mate-1, $n2 mate-2)\n";
}

for my $i (0 .. $n1 - 1)
{
    push(@pe_read_files, [$opt->pe_read_1->[$i], $opt->pe_read_2->[$i]]);
    push(@inputs, $opt->pe_read_1->[$i], $opt->pe_read_2->[$i]);
}

#
# Set output directories
#

my $base = shift;
my $out_dir = shift;

my $int_dir;			# intermediate files

-d $out_dir or mkdir($out_dir) or die "Cannot create $out_dir: $!";

if ($opt->keep_intermediates)
{
    $int_dir = $out_dir;
}
else
{
    $int_dir = File::Temp->newdir(CLEANUP => 1);
}

$base =~ m,/, and die "Output base may not have slash characters\n";

my $output_name = $opt->output_name || $base;
$output_name =~ s/\s+/_/g;

#
# Mapping
# 

# mapping_mode=$(
#     head -n 400 "${input_fastq}" \
#     | awk \
#     -v "thresh=${length_threshold}" \
#     'NR % 4 == 2 {s+= length}END {if (s/(NR/4) > thresh) {print "map-ont"} else {print "sr"}}'
#    )

#
# Determine mapping mode (short or long) based on average read size of the start of the input
#

open(F, "-|", "gunzip", "-c", "-q",  "-d", "-f", $inputs[0]) or die "Cannot open $inputs[0]: $!";

my($count, $total);
while (<F>)
{
    last if $. > 400;
    if ($. % 4 == 2)
    {
	chomp;
	$count++;
	$total += length($_);
    }
}
close(F);
my $avg = $total / $count;

my $mapping_mode = "sr";

if ($avg > $opt->length_threshold || $opt->nanopore)
{
    $mapping_mode = "map-ont";
}

my @minimap_opts = (-K => "20M", 	# Minibatch size
		    '-a',		# Output SAM format alignment
		    -x => $mapping_mode,
		    -t => $opt->threads);

# Trim polyA tail for alignment (33 bases)
my $trimmed = "$int_dir/reference_trimmed.fa";

my $ok = $runner->run(["seqtk", "trimfq", "-e", 33, $reference], '>',  $trimmed);

$ok or die "Failure $? running seqtk\n";

#
# Run the mapper
# 

$runner->run(["minimap2",
	  @minimap_opts,
	  $trimmed,
	  @inputs,
	  "-o", "$int_dir/minimap.out"]);

if ($opt->delete_reads)
{
    print STDERR "Deleting inputs @inputs\n";
    unlink(@inputs);
}

$runner->run_with_timeout(960, ["samtools",
			    "view",
			    "-u",
			    "-h",
			    "-q", $opt->min_quality,
			    "-F", 4,
			    "$int_dir/minimap.out"],
		      '|',
		      ["samtools",
		       "sort",
		       "--threads", $opt->threads,
		       "-o", "$int_dir/$base.sorted.bam",
		       "-"]
		     );
unlink("$int_dir/minimap.out");
$runner->run(["samtools", "index", "$int_dir/$base.sorted.bam"]);

my $ivar_file = "$int_dir/$base.ivar";

$runner->run(["ivar",
	  "trim",
	  "-e",
	  "-q", 0,
	  "-i", "$int_dir/$base.sorted.bam",
	  "-b", "$bed_tmp",
	  "-p", $ivar_file],
	    '>', "$out_dir/$base.primer-trim.txt");

#
# Read the primer trimming report and extract primer-trim.tbl
#
if (open(T, "<", "$out_dir/$base.primer-trim.txt"))
{
    while (<T>)
    {
	last if (/^Primer Name.*Read Count/);
	if (/Found\s+(\d+)\s+mapped/)
	{
	    push(@stats, [mapped_reads => $1]);
	}
	if (/Found\s+(\d+)\s+unmapped/)
	{
	    push(@stats, [unmapped_reads => $1]);
	}
	if (/Found\s+(\d+)\s+primers/)
	{
	    push(@stats, [primer_count => $1]);
	}
    }
    if ($_)
    {
	if (open(O, ">", "$out_dir/$base.primer-trim.tbl"))
	{
	    print O $_;
	    while (<T>)
	    {
		last unless /\t/;
		print O $_;
	    }
	    while (<T>)
	    {
		if (/Trimmed\s+primers\s+from\s+(\S+)%\s+\((\d+)/)
		{
		    push(@stats,
			 [primer_trim_count => $2],
			 [primer_trim_pct => $1]);
		}
	    }
			
	    close(O);
	}
    }
    close(T);
}
    

#
# We are hitting seemingly random timeouts on Bebop with this sort.
#
# Work around it with a timeout and rerun attempt.
#

$runner->run_with_timeout($opt->samtools_sort_timeout, ["samtools",
						    "sort",
						    "$ivar_file.bam",
						    "--threads", $opt->threads,
						    "-o", "$int_dir/$base.isorted.bam"]);

$runner->run(["samtools", "index", "$int_dir/$base.isorted.bam"]);

$runner->run(["samtools",
	  "mpileup",
	  "--fasta-ref", $reference,
	  "--max-depth", $opt->max_depth,
	  "--count-orphans",
	  "--no-BAQ",
	  "--min-BQ", 0,
	  "$int_dir/${base}.isorted.bam"],
	 '>',
	 "$int_dir/$base.pileup");

my $pileup_compress_handle = start(["gzip", "-c", "$int_dir/$base.pileup"],
				   '>',
				   "$out_dir/$base.pileup.gz");

$runner->run(["ivar",
	  "variants",
	  "-p", $ivar_file,
	  "-r", $reference,
	  "-g", reference_gff_path,
	  "-t", 0.6],
	 "<",
	 "$int_dir/$base.pileup");

$runner->run(["ivar",
	  "consensus",
	  "-p", $ivar_file,
	  "-m", $opt->min_depth,
	  "-t", 0.6,
	  "-n", "N"],
	 "<",
	 "$int_dir/$base.pileup");

$runner->run(["sed",
	  '/>/ s/$/ | One Codex consensus sequence/'],
	 "<",
	 "$ivar_file.fa",
	 '|',
	 [qw(seqtk seq -l 60 -)],
	 ">",
	 "$out_dir/$base.fasta");

$runner->run(["cat", $reference, "$out_dir/$base.fasta"],
	 '|',
	 ["mafft", "--auto", "-"],
	 '>',
	 "$out_dir/$base.align");


#
# Create coverage plot
#

$runner->run(["samtools", "depth", "$int_dir/$base.isorted.bam"], '>', "$out_dir/$base.depth");

if (-s "$out_dir/$base.depth")
{
    eval {
	$ENV{GDFONTPATH} = "/usr/share/fonts/liberation";
	my $plot = <<END;
set term png font "LiberationSans-Regular"
set xlabel "Position"
set ylabel "Depth"
set title "Coverage depth for $base"
set output "$out_dir/$base.png"
plot "$out_dir/$base.depth" using 2:3 with impulses title ""
set output
END
    $runner->run(["gnuplot"], "<", \$plot);
    };

    eval {
	$ENV{GDFONTPATH} = "/usr/share/fonts/liberation";
	my $plot = <<END;
set term png font "LiberationSans-Regular"
set yrange [0:250]
set xlabel "Position"
set ylabel "Depth"
set title "Coverage depth for $base"
set output "$out_dir/$base.detail.png"
plot "$out_dir/$base.depth" using 2:3 with impulses title ""
set output
END
    $runner->run(["gnuplot"], "<", \$plot);
    };


}

system("mv", "$int_dir/$base.isorted.bam", "$out_dir/$base.sorted.bam");
system("mv", "$int_dir/$base.isorted.bam.bai", "$out_dir/$base.sorted.bam.bai");
#system("gzip", "-f", "$out_dir/$base.pileup");
system("mv", "$ivar_file.tsv", "$out_dir/$base.variants.tsv");

#
# Compute some statistics
#

open(S, ">", "$out_dir/$base.statistics.tsv") or die "Cannot write $out_dir/$base.statistics.tsv: $!";

eval {
    my($depth_vals) = rcols("$out_dir/$base.depth", 2);

    printf S "depth_mean\t%.1f\n", $depth_vals->avg();
    printf S "depth_median\t%.1f\n", $depth_vals->median();
    printf S "depth_stdv\t%.1f\n", $depth_vals->stdv();
    printf S "depth_min\t%d\n", $depth_vals->min();
    printf S "depth_max\t%d\n", $depth_vals->max();

    open(F, "<", "$out_dir/$base.fasta") or die "Cannot open $out_dir/$base.fasta: $!";
    my($id, $def, $seq) = read_next_fasta_seq(\*F);
    close(F);
    my $nblocks = 0;
    my $ncount = 0;
    while ($seq =~ /([nN]+)/g)
    {
	$nblocks++;
	$ncount += length($1);
    }

    print S "n_count\t$ncount\n";
    print S "n_blocks\t$nblocks\n";
    print S "fasta_length\t" . length($seq) . "\n";

    print S join("\t", @$_) . "\n" foreach @stats;

    open(V, "<", "$out_dir/$base.variants.tsv") or die "Cannot open $out_dir/$base.variants.tsv: $!";
    my $vc = 0;
    $_ = <V>;
    $vc++ while (<V>);
    close(V);
    print S "variant_count\t$vc\n";
    
    
    close(S);
};
if ($@)
{
    warn "Error processing statistics: $@";
}
    

print STDERR "Waiting for pileup gzip to finish\n";
$pileup_compress_handle->finish();

print STDERR  JSON::XS->new->pretty(1)->canonical(1)->encode($runner->report);
