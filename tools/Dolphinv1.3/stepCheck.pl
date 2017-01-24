#!/usr/bin/env perl

#########################################################################################
#                                       stepCheck.pl
#########################################################################################
# 
#  This program checks the input parameters and files. If there is a problem in any of 
#  the params, this program returns an error. 
#  This program also copy/makes links and unzip the input files into the necessary locations 
#  for downstream analysis.
#
#
#########################################################################################
# AUTHORS:
#
# Alper Kucukural, PhD 
# Jul 4, 2014
#########################################################################################


############## LIBRARIES AND PRAGMAS ################

 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 
#################### VARIABLES ######################
 my $input            = "";
 my $resume           = "";
 my $outdir           = "";
 my $jobsubmit        = "";
 my $barcode          = "";
 my $adapter          = ""; 
 my $trim             = "";
 my $downloadS3       = ""; 
 my $username         = "";
 my $config           = "";
 my $wkey             = "";
 my $runparamsid      = "";
 my $servicename      = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $command=$0." ".join(" ",@ARGV); ####command line copy
my $execdir= dirname($0);

GetOptions( 
    'input=s'        => \$input,
    'resume=s'       => \$resume,
    'outdir=s'       => \$outdir,
    'barcode=s'      => \$barcode,
    'username=s'     => \$username,
    'adapter=s'      => \$adapter,
    'trim=s'         => \$trim,
    'wkey=s'         => \$wkey,
    'config=s'       => \$config,
    'downloadS3=s'   => \$downloadS3,
    'paramsid=s'     => \$runparamsid,
    'servicename=s'  => \$servicename,
    'jobsubmit=s'    => \$jobsubmit,
    'help'           => \$help, 
    'version'        => \$print_version,
) or die("Unrecognized options.\nFor help, run this script with -help option.\n");

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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($input eq "") or ($outdir eq "") );	

################### MAIN PROGRAM ####################
#    maps the reads to the ribosome and put the files under $outdir/after_ribosome directory
if (lc($resume) eq "fresh")
{
   `rm -rf $outdir/tmp/track/*`;
   `rm -rf $outdir/seqmapping`;
   `rm -rf $outdir/tophat`;
   `rm -rf $outdir/rsem`;
   `rm -rf $outdir/*tdf*`;
   `rm -rf $outdir/*ucsc*`;
   `rm -rf $outdir/agg`;
   `rm -rf $outdir/macs`;
}

$input=~s/s3:\/\//s3__\/\//g;
$input=~s/:+/:/g;
$input=~s/\s//g;
$input=~s/:$//g;

my @pfiles=split(/:/,$input);

if(lc($adapter) ne "none")
{
   die "Error 66: please check the adapter:$adapter" unless ($adapter =~/^[ACGT:]+$/);
}
if(lc($trim) ne "none")
{
   my @nts=split(/[,\s\t:]+/,$trim);
   foreach my $ntlen (@nts)
   {
      if (length($ntlen)>=1)
      {
         die "Error 67: please check the trim lengths:$ntlen" unless ($ntlen =~/^[\d]+$/);
      }
   }   
}
my $trackdir   = "$outdir/tmp/track";
$outdir   = "$outdir/input";
`mkdir -p $outdir`;
die "Error 15: Cannot create the directory:".$outdir if ($?);

my %prefiles=();
my %s3com=();
foreach my $line (@pfiles)
{ 
  print $line."\n";

  my @files=split(/[,\s\t]+/,$line);
  my $libname="data";
  
  my $maxindex=2;
  my $offset=1;
  if(lc($barcode) eq "none")
  {
    $libname=$files[0];
    print "Libname:[$libname]\n";
    $maxindex=3;
    $offset=0;
  }
  for (my $i=($maxindex-2); $i<@files; $i++)
  {
      my $file=$files[$i];
      print "[$file]\n";
      my $pairedstr="";
      $pairedstr=".".($i+$offset) if (scalar(@files)>($maxindex-1));
      if ($file =~/s3__\/\//){
         $file =~s/s3__\/\//s3:\/\//;
         my @fparts = split(/\//, $file);
         my $outfile = "$outdir/s3tmp/".$fparts[-1];
         my $coms3 = "$downloadS3 -i $file -o $outfile -u $username -c $config && ";
         print $coms3."\n";
         $s3com{$libname.$pairedstr}.=$coms3;
         $prefiles{$libname.$pairedstr}.=$outfile." ";
      }else{
         $prefiles{$libname.$pairedstr}.=$file." ";
         die "Error 64: please check the file:".$file unless (checkFile($file, "$trackdir/".$servicename.$libname.$pairedstr.".success"));
      }
  }
}
my $i=1;
my $com="";

foreach my $libname (keys %prefiles)
{
  my $filestr=$prefiles{$libname};
  if($libname!~/^$/)
  {
    my $cat="cat";
    $cat = "zcat" if ($filestr=~/\.gz/);
    $com = "$cat $filestr > $outdir/$libname.fastq";
    print "[[$filestr]]\n";
    $com = "mkdir -p $outdir/s3tmp/ && ".$s3com{$libname}.$com if ($filestr=~/s3tmp/);
 
    print "COMMAND: [ $com ]\n\n";
    my $job=$jobsubmit." -n ".$servicename.$libname." -c \"$com\"";
    print $job."\n\n";
    `$job`;
	die "Error 25: Cannot run the job:".$job if ($?);
    $i++;
  }
 }

if(lc($barcode) ne "none")
{  
  my @blines=split(/:/,$barcode);
  foreach my $line (@blines)
  {
    if (length($line)>1)
    {
      my @defs=split(/[,\s\t]+/,$line);
      die "Error 65: please check the barcode definitions:$line" unless (scalar(@defs)>1);
    }
  }
}
 

sub checkFile
{
 my ($file, $success_file) = @_;
 print "$success_file\n";
 return 1 if ($file =~/^s3__\/\//);
 return 1 if (-e $success_file);
 return 1 if (-e $file);
 return 0;
}

__END__


=head1 NAME

stepCheck.pl

=head1 SYNOPSIS  

stepCheck.pl -i input <fastq> 
            -o outdir <output directory> 

stepCheck.pl -help

stepCheck.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -i  input file <fastq format> 

fastq files has to be separated with ":". If it is paired end the paired end files has to ber separated by ","

Ex: For single end;

test1.fastq:test2.fastq:ctrl1.fastq:ctrl2.fastq

for paired end;

test1_R1.fastq,test1_R2.fastq:ctrl1_R1.fastq,ctrl1_R2.fastq

    
=head2 -o outdir <output directory>

the output files will be "$outdir/input" 

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program map the reads to rRNAs and put the rest into other files 

=head1 EXAMPLE


stepCheck.pl -i test1.fastq:test2.fastq:ctrl1.fastq:ctrl2.fastq
            -o ~/out
            -b ~/bowtie_dir/bowtie
            -p "-p 8 -n 2 -l 20 -M 1 -a --strata --best"
            -r ~/bowtie_ind/rRNA

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


