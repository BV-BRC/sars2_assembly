
#
# Comprehensive SARS2 Analysis application
#
# This is a coordination/wrapper application that invokes a series of
# other applications to perform the actual work. 
#

use Bio::KBase::AppService::AppScript;
use Bio::P3::SARS2Assembly::ComprehensiveAnalysis;
use strict;
use Data::Dumper;

my $cga = Bio::P3::SARS2Assembly::ComprehensiveAnalysis->new();

my $script = Bio::KBase::AppService::AppScript->new(sub { $cga->run(@_); }, sub { $cga->preflight(@_); });

my $rc = $script->run(\@ARGV);

exit $rc;
