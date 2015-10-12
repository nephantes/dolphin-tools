#!/usr/bin/env perl

#########################################################################################
#                                       stepFastQC.pl
#########################################################################################
# 
#  This program runs fastQC program
#
#########################################################################################
# AUTHORS:
#
# Alper Kucukural, PhD 
# 
#########################################################################################


############## LIBRARIES AND PRAGMAS ################

 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 
#################### VARIABLES ######################
 my $outdir           = "";
 my $barcode          = "";
 my $prog             = "";
 my $jobsubmit        = "";
 my $servicename      = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";

################### PARAMETER PARSING ####################

my $command=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
    'outdir=s'       => \$outdir,
    'prog=s'         => \$prog,
    'barcode=s'      => \$barcode,
    'jobsubmit=s'    => \$jobsubmit,
    'servicename=s'  => \$servicename,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($prog eq "") or ($outdir eq "") or ($barcode eq "") );	

 
################### MAIN PROGRAM ####################
#    maps the reads to the the genome and put the files under $outdir directory


my $inputdir="";
if ($barcode=~/NONE/g)
{
  $inputdir = "$outdir/input";
}
else
{
  $inputdir = "$outdir/seqmapping/barcode";
}
print $inputdir."\n";
$outdir   = "$outdir/fastqc";
`mkdir -p $outdir`;
die "Error 15: Cannot create the directory:".$outdir  if ($?);

my $com="";
$com=`ls $inputdir/*.fastq`;
die "Error 64: please check the if you defined the parameters right:" unless ($com !~/No such file or directory/);

print $com;
my @files = split(/[\n\r\s\t,]+/, $com);

foreach my $file (@files)
{ 
 $file=~/.*\/(.*).fastq/;
 my $bname=$1;
 die "Error 64: please check the file:".$file unless (checkFile($file));
 my $dir=$outdir."/".$bname;
 `mkdir -p $dir`;
 die "Error 15: Cannot create the directory:".$dir if ($?);
 $com="$prog ".$file." -o $dir";
 my $job=$jobsubmit." -n ".$servicename."_".$bname." -c \"$com\"";
 print $job."\n";
 `$job`;
 die "Error 25: Cannot run the job:".$job if ($?);
}


sub checkFile
{
 my ($file) = $_[0];
 return 1 if (-e $file);
 return 0;
}

__END__


=head1 NAME

stepFastQC.pl

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






