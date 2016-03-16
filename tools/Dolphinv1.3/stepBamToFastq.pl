#!/usr/bin/env perl

#########################################################################################
#                                      bamToFastq.pl
#########################################################################################
# 
# Converts bam files to Fastq
#
#########################################################################################
# AUTHORS:
#
# Alper Kucukural, PhD
# 13/03/2015
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
 my $type             = "";
 my $paired           = "";
 my $samtools         = "";
 my $cmd              = "";
 my $jobsubmit        = "";
 my $servicename      = "";
 my $print_version    = "";
 my $help             = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $command=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
    'type=s'         => \$type,
    'outdir=s'       => \$outdir,
    'cmd=s'          => \$cmd,
    'samtools=s'     => \$samtools,
    'paired=s'       => \$paired,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($type eq "") or ($outdir eq "")  );	

################### MAIN PROGRAM ####################
# runs the picard program

my $outd  = "$outdir/seqmapping/$type";

`mkdir -p $outd`;

my @files=();
print $type."\n";
if ($type eq "RSEM")
{ 
   my $indir   = "$outdir/rsem";
   @files = <$indir/pipe*/*.genome.sorted.bam>;
}
elsif ($type eq "chip")
{ 
   my $indir   = "$outdir/seqmapping/chip";
   @files = <$indir/*.sorted.bam>;
}
elsif (lc($type) eq "tophat")
{
   my $indir   = "$outdir/tophat";
   print $indir."\n";
   @files = <$indir/pipe*/*.sorted.bam>;
}
else
{
   my $indir   = "$outdir/$type";
   print $indir."\n";
   @files = <$indir/*.bam>;
}

foreach my $d (@files){ 
  my $dirname=dirname($d);
  my $libname=basename($d, ".bam");
  $libname=~s/\.sorted//g;
 
  my $com= " samtools sort -n $d $outd/$libname.sorted && ";
  $com .= " $cmd -i $outd/$libname.sorted.bam ";

  if (lc($paired) eq "paired" ) {
         $com .= " -fq $outd/$libname.1.fastq ";
         $com .= " -fq2 $outd/$libname.2.fastq && ";
  }
  else
  {
     $com .= " -fq $outd/$libname.fastq && ";
  }
  $com .= " rm -rf $outd/$libname*.bam";
  
  my $job=$jobsubmit." -n ".$servicename."_".$libname." -c \"$com\"";
  print "\n\n=====\n$job\n======\n\n";
  `$job`;
  die "Error 25: Cannot run the job:".$job if ($?);
}

__END__


=head1 NAME

bam2Fastq.pl

=head1 SYNOPSIS  

bamToFastq.pl 
            -o outdir <outdir> 
            -t type <type dedup|merged|tophat>
            -c cmd <bedtols bamToFastq>
            -p paired <paired|no> 

bamToFastq.pl -help

bamToFastq.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -f filename <output directory>

bowtie std output

=head2 -n name <picard running line> 

sample name

=head2  -p <paired>  

paired end library or not <yes|no>

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

This program runs the cufflinks after tophat mappings

=head1 EXAMPLE

bamToFastq.pl 
            -o outdir <outdir> 
            -t type <type dedup|merged|tophat>
            -c cmd <bedtols bamToFastq
            -p paired <paired|no> 

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

