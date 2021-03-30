#
# The SARS2 genome assembly application.
#

=head1 NAME

App-SARS2Assembly - assemble a set of reads

=head1 SYNOPSIS

    App-SARS2Assembly [--preflight] service-url app-definition parameters

=head1 DESCRIPTION

Assemble a set of SARS2 reads.

=head2 PREFLIGHT INFORMATION

On a preflight request, we will generate a JSON object with the following
key/value pairs:

=over 4

=item ram

Requested memory. For standard run we request 128GB.

=item cpu

Requested CPUs.
    
=cut

use strict;
use Cwd;
use Carp;
use Data::Dumper;
use File::Path 'make_path';
use File::Temp;
use File::Basename;
use IO::Handle;
use IPC::Run 'run';
use POSIX;
use File::Slurp;
use JSON::XS;
use gjoseqlib;

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::ReadSet;

our %default_platform_recipe = (illumina => 'cdc-illumina',
				nanopore => 'cdc-nanopore',
				);
#our %default_platform_recipe = (illumina => 'onecodex',
#				nanopore => 'onecodex',
#				);
our %valid_platform_recipe =
    (
     illumina => {
	 'cdc-illumina' => 1,
	 'onecodex' => 1,
     },
     nanopore => {
	 'cdc-nanopore' => 1,
	 'artic-nanopore' => 1,
	 'onecodex' => 1,
     },
    );

#
# requirements are [cpu-count, memory]
#
our %run_requirements = ('cdc-illumina' => [ 12, '16G'],
			 'cdc-nanopore' => [ 8, '16G'],
			 'artic-nanopore' => [ 8, '16G'],
			 'onecodex' => [8, '16G'],
			 );

our %recipe_tool = ('cdc-illumina' => 'sars2-cdc-illumina',
		    'cdc-nanopore' => 'sars2-cdc-nanopore',
		    'artic-nanopore' => 'sars2-artic-nanopore',
		    'onecodex' => 'sars2-onecodex');
#
# Some tools require specific options for some library types.
#
our %tool_library_options;
$tool_library_options{'sars2-onecodex', 'nanopore'} = ['--nanopore'];


my $script = Bio::KBase::AppService::AppScript->new(\&assemble, \&preflight);

my $download_path;

my $rc = $script->run(\@ARGV);

exit $rc;

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    print STDERR "preflight genome ", Dumper($params, $app);

    my $token = $app->token();
    my $ws = $app->workspace();

    my $readset;
    eval {
	$readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params);
    };
    if ($@)
    {
	die "Error parsing assembly parameters: $@";
    }

    my($ok, $errs, $comp_size, $uncomp_size) = $readset->validate($ws);

    if (!$ok)
    {
	die "Reads as defined in parameters failed to validate. Errors:\n\t" . join("\n\t", @$errs);
    }

    my($recipe, $details) = determine_recipe($app_def, $ws, $params, $readset);

    print STDERR "comp=$comp_size uncomp=$uncomp_size recipe=$recipe\n";

    my($est_cpu, $est_ram) = @{$run_requirements{$recipe}};

    my $est_comp = $comp_size + 0.75 * $uncomp_size;
    $est_comp /= 1e6;
    my $est_storage = int(1.3e6 * $est_comp / 0.75);

    my $est_time = 3600 * 3;
    if ($est_comp > 1000)
    {
	$est_time = 86400;
    }

    return {
	cpu => $est_cpu,
	memory => $est_ram,
	runtime => $est_time,
	storage => $est_storage,
    };
}

sub determine_recipe
{
    my($app_def, $ws, $params, $readset) = @_;

    my @libs = grep { ! $_->{derived_from} } $readset->libraries;

    #
    # Count the number of libraries by platform.
    #
    my(%lib_pe_platform, %lib_se_platform);

    my @pe_libs = grep { $_->is_paired_end() } @libs;
    my @se_libs = grep { $_->is_single_end() } @libs;
    
    $lib_pe_platform{$_->{platform}}++ foreach @pe_libs;
    $lib_se_platform{$_->{platform}}++ foreach @se_libs;

    my @pe_platforms = keys %lib_pe_platform;
    my @se_platforms = keys %lib_se_platform;

    if (@pe_platforms && @se_platforms)
    {
	die "SARS2Assembly does not allow a mix of single-end and paired-end read data\n";
    }

    if (@pe_libs > 1)
    {
	die "SARS2Assembly does not allow more than one paired-end read library\n";
    }

    my $lib = @libs[0];

    my $details = {};

    my($recipe_def) = grep { $_->{id} eq 'recipe' } @{$app_def->{parameters}};
    my $valid_recipes = $recipe_def->{enum};

    my $recipe = $params->{recipe};

    my $ok = grep { $_ eq $recipe } @$valid_recipes;
    $ok or die "Recipe $recipe is not valid (acceptable values are @$valid_recipes)";

    my $lib_type;
    my $platform = $lib->{platform};
    if (!$platform || $platform eq 'infer')
    {
	$platform = $lib->{metadata}->{platform_name};
    }
    
    if ($platform =~ /nanopore/i)
    {
	$lib_type = 'nanopore';
    }
    elsif ($platform =~ /ion([_\s])*torrent/i)
    {
	$lib_type = 'nanopore';
    }
    elsif ($platform =~ /illumina/i)
    {
	$lib_type = 'illumina';
    }
    elsif ($recipe eq 'auto')
    {
	my $is_pe = $lib->is_paired_end();
	if (defined($is_pe) && $is_pe)
	{
	    warn "Defaulting unknown paired end library to illumina\n";
	    $lib_type = 'illumina';
	}
	else
	{
	    $lib_type = guess_lib_type($ws, $lib);
	    if (!$lib_type)
	    {
		die "Unknown platform $platform and auto recipe requested; could not guess library type";
	    }
	}
    }

    if ($recipe eq 'auto')
    {
	$recipe = $default_platform_recipe{$lib_type};
    }
    elsif (!$valid_platform_recipe{$lib_type}->{$recipe})
    {
	die "Recipe $recipe is not valid for platform $platform";
    }

    #
    # Artic-nanopore doesn't support multiple read files
    #
    if ($recipe eq 'artic-nanopore' && @se_libs > 1)
    {
	die "The artic-nanopore recipe does not accept multiple read files.\n";
    }

    $details->{platform} = $platform;
    $details->{library_type}  = $lib_type;
    $details->{recipe} = $recipe;

    return ($recipe, $details);
}

sub guess_lib_type
{
    my($ws, $lib) = @_;

    #
    # We're only guessing for single ended libs.
    #
    if (!$lib->is_single_end())
    {
	warn "Can only guess for SE libs\n";
	return undef;
    }

    my $shock;
    eval {
	my $path = $lib->{read_file};
	
	my $info = $ws->get({objects => [$path], metadata_only => 1});
	my $obj = $info->[0];
	my($meta) = @$obj;
	$shock = $meta->[11];
    };
    if (!$shock)
    {
	warn "Error determining file type: $@";
	return undef;
    }

    my $buf = $ws->shock_read_bytes($shock, 0, 1000000);
    my %comp_map = ("\x1f\x8b" => 'gunzip',
		    "BZ" => 'bunzip2');

    my $comp = $comp_map{substr($buf, 0, 2)};
    my $temp;
    my $fh;
    if ($comp)
    {
	$temp = File::Temp->new();
	$fh  = new IO::Handle;
	run([$comp], '<', \$buf, '>', "$temp");
	close($temp);
	open($fh, "<", $temp);
    }
    else
    {
	open($fh, "<", \$buf);
    }
    my $n_reads = 0;
    my $total_len = 0;
    my $max = 0;
    while (<$fh>)
    {
	my $h1 = $_;
	my $dat = <$fh>;
	my $h2 = <$fh>;
	my $q = <$fh>;
	if (!defined($q))
	{
	    # we hit eof in the middle of a record and that's fine.
	    last;
	}
	if ($h1 !~ /^@/ || $h2 !~ /^\+/)
	{
	    warn "Invalid FQ data\n";
	    return undef;
	}
	$n_reads++;
	my $l = length($dat);
	$max = $l if $max < $l;
	$total_len += $l;
    }
    close($fh);
    my $avg = $total_len / $n_reads;
    print "total=$total_len n_reads=$n_reads max=$max avg=$avg\n";
    if ($max < 600)
    {
	return 'illumina';
    }
    else
    {
	return 'nanopore';
    }
      
}
sub assemble
{
    my($app, $app_def, $raw_params, $params) = @_;

    print "Begin assembly ",
    JSON::XS->new->pretty(1)->canonical->encode([$app_def, $raw_params, $params]);

    my $token = $app->token();
    my $ws = $app->workspace();

    my $cleanup = $params->{debug_level} > 0 ? 0 : 1;

    my $tmpdir = File::Temp->newdir( CLEANUP => $cleanup );
    print STDERR "Debug=$params->{debug_level} cleanup=$cleanup tmpdir=$tmpdir\n";
    $download_path = $tmpdir;

    my $asm_out = "$tmpdir/assembly";
    mkdir($asm_out) or die "cannot mkdir $asm_out: $!";
    my $stage_dir = "$tmpdir/staging";
    mkdir($stage_dir) or die "cannot mkdir $tmpdir/staging: $!";
    
    my $readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params, 1);

    my $metadata_dir = "$stage_dir/sra-metadata";
    make_path($metadata_dir);
    my($ok, $errs) = $readset->validate($ws, "$metadata_dir/");

    if (!$ok)
    {
	die "Reads as defined in parameters failed to validate. Errors:\n\t" . join("\n\t", @$errs);
    }

    $readset->localize_libraries($stage_dir);

    $readset->stage_in($ws);

    my($recipe, $details) = determine_recipe($app_def, $ws, $params, $readset);

    #
    # If we are running under Slurm, pick up our memory and CPU limits.
    #
    my $mem = $ENV{P3_ALLOCATED_MEMORY};
    my $cpu = $ENV{P3_ALLOCATED_CPU};

    my $tool = $recipe_tool{$recipe};
    $tool or die "No tool found for recipe $recipe\n";

    my $opts = $tool_library_options{$tool, $details->{library_type}} // [];

    my @params = @$opts;
    
    push(@params, "-j", $cpu) if $cpu;

    # Onecodex & CDC illumina recipes supports min read depth parameter
    if (($recipe eq 'onecodex' || $recipe eq 'cdc-illumina') && $params->{min_depth} > 0)
    {
	push(@params, "--min-depth", $params->{min_depth});
    }

    if (($recipe eq 'onecodex') && $params->{max_depth} > 0)
    {
	push(@params, "--max-depth", $params->{max_depth});
    }

    if ($params->{keep_intermediates})
    {
	push(@params, "--keep-intermediates");
    }

    for my $lib (grep { not $_->{derives} } $readset->libraries)
    {
	if ($lib->is_single_end())
	{
	    push(@params, "--se-read", $lib->{read_path});
	}
	else
	{
	    push(@params, "--pe-read-1", $lib->{read_path_1});
	    push(@params, "--pe-read-2", $lib->{read_path_2});
	}
    }

    my @cmd = ($tool,
	       @params,
	       $params->{output_file},
	       $asm_out); 

    print "Start assembler: @cmd\n";

    my $asm_ok = run(\@cmd);
    my $asm_rc = $?;

    my $output_folder = $app->result_folder();

    #
    # Read the fasta to get some counts.
    #
    if (open(CTG, "<", "$asm_out/$params->{output_file}.fasta"))
    {
	my $nblocks;
	my $ncount;
	while (my($id, $def, $seq) = read_next_fasta_seq(\*CTG))
	{
	    while ($seq =~ /([nN]+)/g)
	    {
		$nblocks++;
		$ncount += length($1);
	    }
	}
	close(CTG);
	$details->{total_ns} = $ncount;
	$details->{n_blocks} = $nblocks;
    }

    open(DETAILS, ">", "$asm_out/assembly-details.json");
    print DETAILS JSON::XS->new->pretty(1)->encode($details);
    close(DETAILS);

    my $type_map = {
	bam => "bam",
	bai => "bai",
	tbi => "tbi",
	csi => "csi",
	vcf => "vcf",
	hdf => "hdf",
	'vcf.gz' => "vcf",
	html => "html",
	png => "png",
	fasta => "contigs",
	run_details => 'txt',
	txt => 'txt',
	log => 'txt',
	stdout => "txt",
	stderr => "txt",
	jpg => 'jpg',
	json => 'json',
	xml => 'xml',
    };

    my $save_path = $asm_out;

    write_dir($ws, $save_path, $output_folder, $type_map);
    $ws->create({ objects => [["$output_folder/sra-metadata", 'folder']] });
    write_dir($ws, $metadata_dir, "$output_folder/sra-metadata", $type_map);
    
    if (!$asm_ok)
    {
	die "Assembler failed with rc=$asm_rc";
    }
}


sub write_dir
{
    my($ws, $save_path, $output_folder, $type_map) = @_;
    if (opendir(DIR, $save_path))
    {
	while (my $f = readdir(DIR))
	{
	    next if $f =~ /^\./;
	    my $path = "$save_path/$f";
	    my($suffix) = $f =~ /\.([^.]+)$/;
	    my $type = $type_map->{$suffix} // 'txt';
	    
	    if (-f $path)
	    {
		$ws->save_file_to_file($path, {}, "$output_folder/$f", $type, 1, 1);
	    }
	    else
	    {
		$ws->upload_folder($path, "$output_folder", { type_map => $type_map });
	    }
	}
	closedir(DIR);
    }
    else
    {
	warn "Cannot opendir $save_path: $!";
    }
}
