#!/usr/bin/env perl

#########################################################################################
#                                       stepTrimmer.pl
#########################################################################################
# 
#  This program trims the reads in the files. 
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
 my $trim             = "";
 my $outdir           = "";
 my $jobsubmit        = "";
 my $spaired          = "";
 my $previous         = ""; 
 my $cmd              = ""; 
 my $servicename      = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
        'trim=s'         => \$trim,
	'outdir=s'       => \$outdir,
        'dspaired=s'     => \$spaired,
        'previous=s'     => \$previous,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($trim eq "") or ($outdir eq "") );	

################### MAIN PROGRAM ####################
#    maps the reads to the ribosome and put the files under $outdir/after_ribosome directory

my $inputdir="";
print "$previous\n";
if ($previous=~/NONE/g)
{
  $inputdir = "$outdir/input";
}
else
{
  $inputdir = "$outdir/recmapping/".lc($previous);
}

$outdir   = "$outdir/recmapping/trim";
`mkdir -p $outdir`;
my $com="";
if ($spaired eq "single")
{
 $com=`ls $inputdir/*.fastq`;
}
else
{
 $com=`ls $inputdir/*.1.fastq`;
}

print $com;
my @files = split(/[\n\r\s\t,]+/, $com);

foreach my $file (@files)
{
 die "Error 64: please check the file:".$file unless (checkFile($file));
 if ($spaired eq "single")
 {
    $file=~/.*\/(.*).fastq/;
    my $bname=$1;
    trimFiles($file, $trim, $bname);
 }
 else
 {
    print "PAIRED\n\n";
    $file=~/(.*\/(.*)).1.fastq/;
    my $bname=$2;
    my $file2=$1.".2.fastq";
    die "Error 64: please check the file:".$file2 unless (checkFile($file2));
    $trim=~/([\d]*,[\d]*),([\d]*,[\d]*)/;
    my $trim1=$1;
    my $trim2=$2;
    trimFiles($file, $trim1, $bname.".1");
    trimFiles($file2, $trim2, $bname.".2");
 }
}

sub trimFiles
{
  my ($file, $trim, $bname)=@_;
    my @nts=split(/[,\s\t]+/,$trim);
    print "\n$bname\n\n";
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="-f";
    foreach my $nt (@nts)
    {
      if ($nt!~/^$/ && $nt>0)
      {
        $nt++;
        if ($i%2==0)
        {
           $param="-t";
           $nt--;
        }

        $outfile="$outdir/$bname.fastq.$i.tmp";
        $com.="$cmd -v $param $nt -o $outfile -i $file;";
        $file=$outfile;
      }
      $i++;
    }
    $com.="mv $outfile $outdir/$bname.fastq;rm -rf $outdir/*.tmp";
    print $com."\n";
    `$com`;
    my $job=$jobsubmit." -n ".$servicename."_".$bname." -c \"$com\"";
    #print $job."\n";   
    #`$job`;
}


sub checkFile
{
 my ($file) = $_[0];
 return 1 if (-e $file);
 return 0;
}

__END__


=head1 NAME

stepTrimmer.pl

=head1 SYNOPSIS  

stepTrimmer.pl -i input <fastq> 
            -o outdir <output directory> 
            -b bowtieCmd <bowtie dir and file> 
            -p params <bowtie params> 
            -r ribosomeInd <ribosome Index file>

stepTrimmer.pl -help

stepTrimmer.pl -version

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


stepTrimmer.pl -i test1.fastq:test2.fastq:ctrl1.fastq:ctrl2.fastq
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


