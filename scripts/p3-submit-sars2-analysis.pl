=head1 Submit a PATRIC SARS-CoV-2 Analysis Job

    p3-submit-sars2-analysis [options] output-path output-name

Submit a set of one or more read libraries to the PATRIC SARS-CoV-2 assembly service.

=head1 Usage synopsis

    p3-submit-sars2-analysis [-h] output-path output-name

	Submit an assembly job with output written to output-path and named
	output-name.    

    	The following options describe the inputs to the assembly:

           --workspace-path-prefix STR   Prefix for workspace pathnames as given
    	   			   	 to library parameters.
           --workspace-upload-path STR	 If local pathnames are given as library parameters,
    					 upload the
    					 files to this directory in the workspace.
           --overwrite			 If a file to be uploaded already exists in
    	   				 the workspace, overwrite it on upload. Otherwise
    					 we will not continue the service submission.
    	   --paired-end-lib P1 P2	 A paired end read library. May be repeated.
           --single-end-lib LIB	 	 A single end read library. May be repeated.
	   --srr-id STR		 	 Sequence Read Archive Run ID. May be repeated.

    	The following options describe the processing requested:

       	   --scientific-name STR 	 Scientific name for this genome.
           --taxonomy-id STR	         Taxonomy id for this genome.
           --recipe                      Assembly recipe. Defaults to auto.
	   				 Valid values are cdc-nanopore, cdc-illumina, or artic-nanopore
           --container-id		 (Internal) Use the specified container to run this application

    The following options describe the read libraries:
	
           --platform STR		 The sequencing platform for the next read 
    					 library or libraries. Valid values are
					 infer, illumina, pacbio, nanopore, iontorrent

=cut

use strict;
use Getopt::Long;
use Bio::P3::Workspace::WorkspaceClientExt;
use Bio::KBase::AppService::Client;
use P3AuthToken;
use Try::Tiny;
use IO::Handle;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

use File::Basename;
use JSON::XS;
use Pod::Usage;
use Fcntl ':mode';

#
# We start with a bare parameters struct and fill
# it in as we process the command line parameters.
#

my $params = {
    paired_end_libs => [],
    single_end_libs => [],
    srr_ids => [],
    recipe => "auto",
    scientific_name => "Severe acute respiratory syndrome coronavirus 2",
    taxonomy_id => 2697049,
    domain => "Viruses",
    code => 1,
    input_type => 'reads',
};

my @cur_pair = ();
my $platform = "";

#
# Perhaps default this to the user's home.
#
my $workspace_path_prefix;
my $workspace_upload_path;
my $overwrite;
my $dry_run;

my $token = P3AuthToken->new();
if (!$token->token())
{
    die "You must be logged in to PATRIC via the p3-login command to submit assembly jobs\n";
}
my $ws = Bio::P3::Workspace::WorkspaceClientExt->new();
my $app_service = Bio::KBase::AppService::Client->new();

GetOptions("recipe=s" => \$params->{recipe},
	   "platform=s" => \$platform,
	   "container-id=s" => \$params->{container_id},
           'paired-end-lib=s{2}',  sub {
	       my ($opt_name, $opt_value) = @_;

	       #
	       # Check if the option value looks like another parameter
	       # (starts with --). If so we probably missed a file in the pair.
	       if (@cur_pair == 2)
	       {
		   die "Invalid pair specified\n";
	       }
	       if ($opt_value =~ /^-/)
	       {
		   die "Invalid pair specified with $opt_value\n";
	       }
	       push(@cur_pair, process_filename($opt_value));
	       if (@cur_pair == 2)
	       {
		   my $lib = {
		       read1 => $cur_pair[0],
		       read2 => $cur_pair[1],
		       platform => $platform,
		   };
		   push(@{$params->{paired_end_libs}}, $lib);
		   @cur_pair = ();
	       }
	   },
	   "single-end-lib=s" => sub {
	       my($opt_name, $opt_value) = @_;
		   my $lib = {
		       read => process_filename($opt_value),
		       platform => $platform,
		   };
		   push(@{$params->{single_end_libs}}, $lib);
	   },
	   "srr-id=s" => sub {
	       my($opt_name, $opt_value) = @_;
	       push(@{$params->{srr_ids}}, $opt_value);
	   },
	   "overwrite" => \$overwrite,
	   "dry-run" => \$dry_run,
	   "workspace-path-prefix=s" => \$workspace_path_prefix,
	   "workspace-upload-path=s" => sub {
	       my($opt_name, $opt_value) = @_;
	       #
	       # Ensure the upload path exists and is a directory
	       #
	       my $stat = $ws->stat($opt_value);
	       if (!$stat)
	       {
		   die "Workspace upload path $opt_value does not exist\n";
	       }
	       if (!S_ISDIR($stat->mode))
	       {
		   die "Workspace uplaod path $opt_value is not a folder\n";
	       }
	       $workspace_upload_path = $opt_value;
	   },
	   "help|h" => sub {
	       print pod2usage(-sections => 'Usage synopsis', -verbose => 99, -exitval => 0);
	   },
	   );

@ARGV == 2 or pod2usage(-sections => 'Usage synopsis', -verbose => 99, -exitval => 1);

my $output_path = shift;
my $output_name = shift;

#
# we assume output is in workspace, so just clip the prefix.
$output_path =~ s/^ws://;

if ($output_path !~ m,^/,,)
{
    $output_path = $workspace_path_prefix . "/" . $output_path;
}

#
# Make sure it is a folder.
#
my $stat = $ws->stat($output_path);
if (!$stat || !S_ISDIR($stat->mode))
{
    die "Output path $output_path does not exist\n";
}

$params->{output_path} = $output_path;
$params->{output_file} = $output_name;
my @upload_queue;

sub process_filename
{
    my($path) = @_;
    my $wspath;
    if ($path =~ /^ws:(.*)/)
    {
	$wspath = $1;
	if ($wspath !~ m,^/,)
	{
	    if (!$workspace_path_prefix)
	    {
		die "Cannot process $path: no workspace path prefix set (--workspace-path-prefix parameter)\n";
	    }
	    $wspath = $workspace_path_prefix . "/" . $wspath;
	}
	my $stat = $ws->stat($wspath);
	if (!$stat || !S_ISREG($stat->mode))
	{
	    die "Workspace path $wspath not found\n";
	}
    }
    else
    {
	if (!-f $path)
	{
	    die "Local file $path does not exist\n";
	}
	if (!$workspace_upload_path)
	{
	    die "Upload was requested for $path but an upload path was not specified via --workspace-upload-path\n";
	}
	my $file = basename($path);
	$wspath = $workspace_upload_path . "/" . $file;

	if (!$overwrite && $ws->stat($wspath))
	{
	    die "Target path $wspath already exists and --overwrite not specified\n";
	}

	push(@upload_queue, [$path, $wspath]);
    }
    return $wspath;
}

my $errors;
for my $ent (@upload_queue)
{
    my($path, $wspath) = @$ent;
    my $size = -s $path;
    printf "Uploading $path to $wspath (%s)...\n", format_size($size);
    my $res;
    eval {
	$res = $ws->save_file_to_file($path, {}, $wspath, 'reads', $overwrite, 1, $token->token());
    };
    if ($@)
    {
	die "Failure uploading $path to $wspath\n";
    }
    my $stat = $ws->stat($wspath);
    if (!$stat)
    {
	print "Error uploading (file was not present after upload)\n";
	$errors++;
    }
    elsif ($stat->size != $size)
    {
	printf "Error uploading (filesize at workspace (%s) did not match original size $size)\n",
		$stat->size;
	$errors++;
    }
    else
    {
	print "done\n";
    }
}

die "Exiting due to errors\n" if $errors;

if ($dry_run)
{
    print "Would submit with data:\n";
    print JSON::XS->new->pretty(1)->encode($params);
}
else
{
    my $app = "ComprehensiveSARS2Analysis";
    my $task = eval { $app_service->start_app($app, $params, '') };
    if ($@)
    {
	die "Error submitting analysis to service:\n$@\n";
    }
    elsif (!$task)
    {
	die "Error submitting analysis to service (unknown error)\n";
    }
    print "Submitted analysis with id $task->{id}\n";
}

sub format_size
{
    my($s) = @_;
    return sprintf "%.1f Gbytes", $s / 1e9 if ($s > 1e9);
    return sprintf "%.1f Mbytes", $s / 1e6 if ($s > 1e6);
    return sprintf "%.1f Kbytes", $s / 1e3 if ($s > 1e3);
    return "$s bytes";
}
