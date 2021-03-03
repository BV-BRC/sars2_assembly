package Bio::P3::CmdRunner;

#
# Little module to wrap calls to IPC::Run and to
# log execution and runtimes. We also collect the
# set of commands run, map to full paths, and collect
# version info when desired.
#

use strict;
use Time::HiRes 'gettimeofday';
use IPC::Run;
use Data::Dumper;
use Carp::Always;
use File::Which;
use File::Basename;
use Digest::MD5 'md5_hex';
use Cwd;

sub new
{
    my($class) = @_;
    my $self = {
	commands => {},
	log => [],
    };

    return bless $self, $class;
}

sub make_command_log
{
    my($self, $cmds) = @_;
    my $logged = [];
    for my $c (@$cmds)
    {
	if (ref($c) eq 'ARRAY')
	{
	    $self->{commands}->{$c->[0]}++;
	    push(@$logged, [map { ref($_) ? "$_" : $_ } @$c]);
	    print STDERR "\t@$c\n";
	}
	elsif (ref($c) eq 'SCALAR')
	{
	    my $x = "<scalar ref of length " . length($$c) . ">";
	    push(@$logged, $x);
	    print STDERR "\t$x\n";
	}
	elsif (ref($c))
	{
	    push(@$logged, "$c");
	    print STDERR "\t$c\n";
	}
	else
	{
	    push(@$logged, $c);
	    print STDERR "\t$c\n";
	}
    }
    return $logged;
}

sub run
{
    my($self, @cmds) = @_;

    my $logged = $self->make_command_log(\@cmds);

    my $start = gettimeofday;
    my $ok = IPC::Run::run(@cmds);
    my $end = gettimeofday;
    my $elap = $end - $start;
    print STDERR "Run returns $ok $? elapsed=$elap\n";

    push(@{$self->{log}}, [$logged, getcwd, $start, $end, $elap]);
    $ok or die "Failed running pipeline: \n" . Dumper(\@cmds);
    
}


sub run_with_timeout
{
    my($self, $timeout, @cmds) = @_;

    print STDERR "Execute with timeout=$timeout:\n";
    my $logged = $self->make_command_log(\@cmds);

    while (1)
    {
	my $start = gettimeofday;
	my $ok = eval { IPC::Run::run(@cmds, IPC::Run::timeout($timeout)); };
	my $end = gettimeofday;
	my $elap = $end - $start;
	print STDERR "Run returns $ok $? elapsed=$elap\n";
	if ($@ =~ /IPC::Run/)
	{
	    warn "Run failed with IPC::Run error (retrying): $@";
	    next;
	}
	if (! $ok)
	{
	    die "Failed running pipeline: \n" . Dumper(\@cmds);
	}
	else
	{
	    push(@{$self->{log}}, [$logged, getcwd, $start, $end, $elap]);
	    last;
	}
    }
}

sub get_version
{
    my($self, $cmd) = @_;

    my $base = basename($cmd);
    if ($base eq 'ivar')
    {
	my $v;
	IPC::Run::run([$cmd, "version"], '>', \$v, '<', '/dev/null');
	if ($v =~ /version\s+(\S+)/)
	{
	    return $1;
	}
	return undef;
    }
    elsif ($base eq 'seqtk' ||
	   $base eq 'samtools' ||
	   $base eq 'tabix' ||
	   $base eq 'bgzip' ||
	   $base eq 'bcftools' ||
	   $base eq 'bamclipper.sh')
    {
	my $v;
	IPC::Run::run([$cmd], '2>', \$v, '<', '/dev/null');
	if ($v =~ /Version:\s+(\S+)/)
	{
	    return $1;
	}
	return undef;
    }
    elsif ($base eq 'cutadapt')
    {
	my $v;
	IPC::Run::run([$cmd, "--version"], '>', \$v, '<', '/dev/null');
	if ($v =~ /^(\S+)/)
	{
	    return $1;
	}
	return undef;
    }
    elsif ($base eq 'mafft')
    {
	my $v;
	IPC::Run::run([$cmd, "--version"], '2>', \$v, '<', '/dev/null');
	if ($v =~ /^(\S+)/)
	{
	    return $1;
	}
	return undef;
    }
    elsif ($base eq 'minimap2' ||
	   $base eq 'gnuplot'
	  )
    {
	my $v;
	IPC::Run::run([$cmd, '--version'], '>', \$v, '<', '/dev/null');
	chomp $v;
	return $v;
    }
    elsif ($base eq 'cat' ||
	   $base eq 'sed' ||
	   $base eq 'tee' ||
	   $base eq 'medaka' ||
	   $base eq 'bowtie2-build' ||
	   $base eq 'bowtie2' ||
	   $base eq 'vcf_mask_lowcoverage.pl')
    {
	my $v;
	IPC::Run::run([$cmd, '--version'], '>', \$v, '<', '/dev/null');
	if ($v =~ /(\S+)$/m)
	{
	    return $1;
	}
	return undef;

    }

    return undef;
}

sub report
{
    my($self) = @_;
    my $report = {
	commands => [],
	log => $self->{log},
    };
    for my $cmd (sort keys %{$self->{commands}})
    {
	my $n = $self->{commands}->{$cmd};
	my $v = $self->get_version($cmd);
	my $path = which($cmd);

	my $md = Digest::MD5->new;
	if (open(my $fh, "<", $path))
	{
	    $md->addfile($fh);
	}
	else
	{
	    warn "Cannot open $cmd: $!";
	};
	
	push(@{$report->{commands}}, {
	    command => $cmd,
	    path => $path,
	    version => $v,
	    run_count => $n,
	    md5sum => $md->hexdigest,
	});
    }
    return $report;
}

1;
