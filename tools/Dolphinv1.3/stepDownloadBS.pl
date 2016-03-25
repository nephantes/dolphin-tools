#!/usr/bin/env perl
##########################################################################################
#                                       stepDownloadBS.pl
#########################################################################################
# 
#  This program  downloads files from basespace with R.
#
#
#########################################################################################
# AUTHORS:
#
# Alper Kucukural, PhD 
# Mar 24, 2016
#########################################################################################


############## LIBRARIES AND PRAGMAS ################

 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 

#################### VARIABLES ######################
 my $token            = "";
 my $outdir           = "";
 my $projectID        = "";
 my $pattern          = "";
 my $rscriptCMD       = "Rscript"; 
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
    'token=s'        => \$token,
    'outdir=s'       => \$outdir,
    'id=s'           => \$projectID,
    'pattern=s'      => \$pattern,
    'rscriptCMD=s'   => \$rscriptCMD,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($token eq "") or ($projectID eq "") or ($pattern eq ""));	

################### MAIN PROGRAM ####################
#    downloads the files from basespace
`mkdir -p $outdir`;
downloadFiles($token, $projectID, $pattern, $outdir);

sub downloadFiles
{
my ($token, $projectID, $pattern, $outdir)=@_;
my $output = "$outdir/rscript_downloadBS.R";
my $sessioninfo = "$outdir/sessionInfo.txt";

open(OUT, ">$output");
my $rscript = qq/
library(BaseSpaceR)
ACCESS_TOKEN <- "$token"
PROJECT_ID <- "$projectID"  ## Get proj ID from url of the project

aAuth<- AppAuth(access_token = ACCESS_TOKEN)
selProj <- Projects(aAuth, id = PROJECT_ID, simplify = TRUE) 
sampl <- listSamples(selProj, limit= 1000)
inSample <- Samples(aAuth, id = Id(sampl), simplify = TRUE)
for(s in inSample){ 
    f <- listFiles(s,pattern="$pattern", Extensions = ".gz")
    a<-Name(f)
    down<-0
    if (grepl("^$pattern", a[1])){
         down<-1
         print(Name(f))
    }
    if(down){
      getFiles(aAuth, id= Id(f), destDir = "$outdir\/", verbose = TRUE)
    }
}
/;

print $rscript; 

print OUT $rscript; 
close(OUT);

my $com="$rscriptCMD $output > $sessioninfo 2>&1";
`$com`;
}

__END__

=head1 NAME

stepDownloadBS.pl

=head1 SYNOPSIS  

stepDownloadBS.pl -t <token> 
                  -o <output directory>
                  -i <projectID>
                  -p <pattern> 

stepDownloadBS.pl -help

stepDownloadBS.pl -version

For help, run this script with -help option.

=head1 OPTIONS
    
=head2 -o outdir <output directory>

the output files will be "$outdir" 

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program downloads files from basespace

=head1 EXAMPLE


stepDownloadBS.pl -t ACASCASFQQ324254ACAAD
                  -o ~/outdir
                  -i 121231231
                  -p "E91"

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



