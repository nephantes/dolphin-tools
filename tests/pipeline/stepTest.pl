#!/usr/bin/env perl

#########################################################################################
#                                       stepTest.pl
#########################################################################################
# 
#  This program tests the workflow engine
#
#
#########################################################################################
# AUTHORS:
#
# Alper Kucukural, PhD 
# Jan 19, 2016
#########################################################################################


############## LIBRARIES AND PRAGMAS ################

 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 

#################### VARIABLES ######################
 my $cmd              = "";
 my $outdir           = "";
 my $servicename      = "";
 my $num              = "";
 my $jobsubmit        = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "stepTest v1.0.0";

################### PARAMETER PARSING ####################

my $command=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
	'outdir=s'       => \$outdir,
    'cmd=s'          => \$cmd,
	'num=s'          => \$num,
    'servicename=s'  => \$servicename,
    'jobsubmit=s'    => \$jobsubmit,
	'help'           => \$help, 
	'version'        => \$print_version,
) or die("Unrecognized optioins.\nFor help, run this script with -help option.\n");

if($help){
    pod2usage( {
		'-verbose' => 2, 
		'-exitval' => 1,
	} );
}

if($print_version){
  print "Version ".$version."\n";
  exit;
}

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($cmd eq "") or ($outdir eq "") );	

################### MAIN PROGRAM ####################
#  Prepare count and summary files


for (my $i=0; $i<$num; $i++)
{
	
  my $outd   = "$outdir/$servicename.$i";
  `mkdir -p $outd`;
  die "Error 15: Cannot create the directory:".$outd if ($?);
  $cmd =~s/\"/\\\"/g;
  my $job=$jobsubmit ." -n ".$servicename.".".$i." -c \"$cmd>$outd/std.out\"";
  print $job."\n";
  `$job`;
  die "Error 25: Cannot run the job:".$job if ($?);
}
 

__END__


=head1 NAME

stepTest.pl

=head1 SYNOPSIS  

stepTest.pl -c command <command> 
            -o outdir <output directory> 

stepTest.pl -help

stepTest.pl -version

For help, run this script with -help option.

=head1 OPTIONS
    
=head2 -o outdir <output directory>

the output files will be "$outdir/$servicename.$num"

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program runs a test workflow

=head1 AUTHORS

 Alper Kucukural, PhD

 
=head1 LICENSE AND COPYING

 This program is free software; you can redistribute it and / or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, a copy is available at
 http://www.gnu.org/licenses/licenses.html


