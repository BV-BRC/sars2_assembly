#
# Run a chunk of jobs on a bebop node.
#
# We look up the data based on the SLURM_ARRAY_TASK_ID
# and SLURM_ARRAY_TASK_STEP values; that is, we will
# compute on items task_id thru (task_id + task_step - 1)
# in parallel, with a thread count for each of int(36 / task_step).
#
# For each item, we use p3-sra to download the data to /scratch/item.$n
# We split the output directory based on last two digits of SRA number $L2
# 
# Then we run sars2-onecodex $dir/*fastq $job_output/$L2 $sra_number
#
# bebop-run-chunk sra-def-file output-dir
#

use strict;
use Data::Dumper;
use File::Path 'make_path';

@ARGV == 3 or die "Usage: $0 job-offset sra-def-file output-dir\n";

my $job_offset = shift;
my $sra_defs = shift;
my $output = shift;

my $scratch = "/scratch";

$job_offset =~ /^\d+$/ or die "Invalid job offset '$job_offset'\n";
-d $output or die "Output directory $output does not exist\n";

my $start = $job_offset + $ENV{SLURM_ARRAY_TASK_ID};
my $end = $start + $ENV{SLURM_ARRAY_TASK_STEP} - 1;

my $threads = int(36 / $ENV{SLURM_ARRAY_TASK_STEP});
$threads = 4 if $threads > 4;

print STDERR "Run from $start to $end\n";

open(S, "<", $sra_defs) or die "Cannot open SRA defs file $sra_defs: $!";

my @to_run;
while (<S>)
{
    chomp;
    if ($. >= $start && $. <= $end)
    {
	my($id) = split(/\t/);
	push(@to_run, [$., $id]);
    }
}
close(S);

print STDERR "Run " . Dumper(\@to_run);
    
my @pids;

for my $item (@to_run)
{
    my $pid = fork;
    if ($pid == 0)
    {
	run_one(@$item);
	exit 0;
    }
    push(@pids, $pid);
}

for my $pid (@pids)
{
    print STDERR "Wait for $pid\n";
    my $rc = waitpid($pid, 0);
    print  STDERR "Pid $pid is done with $?\n";
}

sub run_one
{
    my($id, $sra) = @_;
    $ENV{TMPDIR} = $scratch;
    print STDERR "Running $id: $sra\n";
    my $tmp = "$scratch/task-$id";
    make_path($tmp);
    my $rc = system("p3-sra", "--id", $sra, "--out", $tmp);
    $rc == 0 or die "p3_sra failed with $rc for $sra\n";
    my $l2 = substr($sra, 0, 7);
    my $out = "$output/$l2/$sra";
    make_path($out);

    my @fastq = glob("$tmp/*.fastq");
    print STDERR "Run on @fastq\n";
    if (@fastq != 1 && @fastq != 2)
    {
	die "Invalid fastq count @fastq\n";
    }
    my @cmd = ("sars2-onecodex", @fastq, $sra, $out, "--threads", $threads);
    print STDERR "Run in pid $$ for id $id $sra: @cmd\n";
    $rc = system(@cmd);
}
