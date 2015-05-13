#!/usr/bin/env perl

#########################################################################################
#                                       stepBarcode.pl
#########################################################################################
# 
#  This program removes barcode sequence. 
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
 my $barcode          = "";
 my $outdir           = "";
 my $jobsubmit        = "";
 my $spaired          = "";
 my $input            = ""; 
 my $cmd              = ""; 
 my $servicename      = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

GetOptions( 
        'input=s'        => \$input,
        'barcode=s'      => \$barcode,
	'outdir=s'       => \$outdir,
        'dspaired=s'     => \$spaired,
        'cmd=s'          => \$cmd,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($barcode eq "") or ($outdir eq "") );	

################### MAIN PROGRAM ####################
#    maps the reads to the ribosome and put the files under $outdir/after_ribosome directory

my $inputdir   = "$outdir/input";

$outdir   = "$outdir/seqmapping/barcode";
`mkdir -p $outdir`;
my @names=();
my @barcodes=();
open(OUT, ">$outdir/barcode.fa");
$barcode=~s/[,\s]+/\t/g;
$barcode=~s/:+/\n/g;
my @nms=split(/\n/,$barcode);
foreach my $n(@nms)
{
  $n=~/(.*)\t(.*)/;
  if ($1!~/^Distance/ || $1!~/^Format/)
  {
    push(@names, $1);
    push(@barcodes, $2);
  }
}
print OUT "$barcode";
close(OUT);

my $com="";
my @files=split(/:/,$input);
my @cmds=split(/:/,$cmd);
my $cmdSE=$cmds[0];
my $cmdPE=$cmds[1];

 my $bname="";
 if ($spaired eq "single")
 {
    my $filename=getFileName($inputdir, 0);

    print $filename."\n\n";
    $filename=~/.*\/(.*)\.(.*)/;
    my $nm=$1;
    my $ext=$2;
    $bname=$nm;
    my $mvcom="";
    for (my $i=0; $i<@names; $i++)
    {
      if ($names[$i]!~/^$/){
       $mvcom.="mv $outdir/$barcodes[$i]/$nm.$ext $outdir/$names[$i].fastq;";
       $mvcom.="rmdir $outdir/$barcodes[$i];";
      }
    }
    $com="$cmdPE -b $outdir/barcode.fa -f $filename -d $outdir > /dev/null;$mvcom";  
 }
 else
 {
    print "PAIRED\n\n";
    
      my $file1=getFileName($inputdir, 1);
      my $file2=getFileName($inputdir, 2);
      print "$file1:$file2\n\n";
      $file1=~/.*\/(.*)\.(.*)/;
      my $nm1=$1;
      my $ext1=$2;
      $bname=$nm1;
      $file2=~/.*\/(.*)\.(.*)/;
      my $nm2=$1;
      my $ext2=$2;
      my $mvcom="";
      for (my $i=0; $i<@names; $i++)
      {
        if ($names[$i]!~/^$/){
         $mvcom.="mv $outdir/$barcodes[$i]/$nm1.$ext1 $outdir/$names[$i].1.fastq;";
         $mvcom.="mv $outdir/$barcodes[$i]/$nm2.$ext2 $outdir/$names[$i].2.fastq;";
         $mvcom.="rmdir $outdir/$barcodes[$i];";
        }
      }

      #$com="$cmdPE -bcfile $outdir/barcode.fa -in $file1 -pair2File $file2 -outdir $outdir > /dev/null;$mvcom";  
      $com="$cmdPE -b $outdir/barcode.fa -f $file1 $file2 -d $outdir > /dev/null;$mvcom";  
      #$com="$mvcom";  
 }
 #print $com."\n\n";
 #`$com`;
 
 my $job=$jobsubmit." -n ".$servicename."_".$bname." -c \"$com\"";
 print $job."\n";   
 `$job`;

sub getFileName
{
  my ($inputdir, $i)= @_;
  
  my  $pairstr="";
  $pairstr=".$i" if ($i>0);
  my $filename="$inputdir/data$pairstr.fastq";
  die "Error 64: please check the file or check if you chose single and paired end library right!".$filename unless (checkFile($filename));
  return $filename;
}

sub checkFile
{
 my ($file) = $_[0];
 return 1 if (-e $file);
 return 0;
}

__END__


=head1 NAME

stepBarcode.pl

=head1 SYNOPSIS  

stepBarcode.pl -i input <fastq> 
            -o outdir <output directory> 
            -b bowtieCmd <bowtie dir and file> 
            -p params <bowtie params> 
            -r ribosomeInd <ribosome Index file>

stepBarcode.pl -help

stepBarcode.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -i  input file <fastq format> 

fastq files has to be separated with ":". If it is paired end the paired end files has to ber separated by ","

Ex: For single end;

test1.fastq:test2.fastq:ctrl1.fastq:ctrl2.fastq

for paired end;

test1_R1.fastq,test1_R2.fastq:ctrl1_R1.fastq,ctrl1_R2.fastq

    
=head2 -o outdir <output directory>

the output files will be "$outdir/after_ribosome" 

=head2 -b bowtieCmd <bowtie dir and file> 

Fullpath of bowtie executable file. Ex: ~/bowtie_dir/bowtie

=head2  -p params <bowtie params> 

Bowtie running parameteres. Ex: "-p 8 -n 2 -l 20 -M 1 -a --strata --best"

=head2  -r ribosomeInd <ribosome Index file>

Ribosomal index files. Ex: ~/bowtie_ind/rRNA


=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program map the reads to rRNAs and put the rest into other files 

=head1 EXAMPLE


stepBarcode.pl -i test1.fastq:test2.fastq:ctrl1.fastq:ctrl2.fastq
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


