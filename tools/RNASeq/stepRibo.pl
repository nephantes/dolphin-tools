#!/usr/bin/env perl

#########################################################################################
#                                       stepRibo.pl
#########################################################################################
# 
#  This program maps the reads to Ribosomal RNAs, if there is no Ribosomal RNA for this
#  genome the program copy the input file into the necessary location for downstream 
#  analysis.
#
#
#########################################################################################
# AUTHORS:
#
# Hennady Shulha, PhD 
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
 my $input            = "";
 my $outdir           = "";
 my $jobsubmit        = "";
 my $ribosomeInd      = "";
 my $bowtiecmd        = ""; 
 my $servicename      = "";
 my $param            = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
	'input=s'        => \$input,
	'outdir=s'       => \$outdir,
        'bowtieCmd=s'    => \$bowtiecmd,
        'ribosomeInd=s'  => \$ribosomeInd,
        'servicename=s'  => \$servicename,
        'jobsubmit=s'    => \$jobsubmit,
        'param=s'        => \$param,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($input eq "") or ($outdir eq "") or ($bowtiecmd eq "") );	

 print "[$servicename]\n";
################### MAIN PROGRAM ####################
#    maps the reads to the ribosome and put the files under $outdir/after_ribosome directory

$outdir   = "$outdir/after_ribosome";

mkdir $outdir if (! -e $outdir);

my @prefiles=split(/:/,$input);
  
for(my $i=0;$i<@prefiles;$i++) 
{
  my @files=split(/,/,$prefiles[$i]);
  if (@files>=3)
  {
    my $libname=$files[1];
    my $str_file=$files[2];
    my $com="$bowtiecmd --un $outdir/$libname.notR -x $ribosomeInd $str_file --al $outdir/$libname.yesR  >  $outdir/$libname.mapped\n";
      
    
    if ((-s $files[2])==0 && $files[2]!~/\*/)
    {
        print "ERROR:Please check the file:".$files[2].", and try again!!!\n";
        exit 1;
   
    }
        
    if ($files[2]=~/\*/) {
	mapMultipleFiles(\@files, $jobsubmit, $servicename, $bowtiecmd, $ribosomeInd, $outdir);
    }
    else
    {
     if (@files==4)
     {
      
      $str_file="-1 ".$files[2]." -2 ".$files[3];
      $com="$bowtiecmd --un-conc $outdir/$libname.notR -x $ribosomeInd $str_file --al-conc $outdir/$libname.yesR  >  $outdir/$libname.mapped\n";
      
      if ((-s $files[2])==0)
      {
        print "ERROR:Please check the file:".$files[3].", and try again!!!\n";
        exit 1;
      }
     }
  
     if ($ribosomeInd!~/^$/)
     {
        print "[".$com."]\n\n";
       print "[".@files.":".@prefiles.":".$jobsubmit."]\n";
       if(((@files-2)+@prefiles)>2 && $jobsubmit!~/^$/)
       {
         my $job=$jobsubmit." -n ".$servicename."_".$libname." -c \"$com\"";
         #print $job."\n";   
         `$job`;
	 print "THERE\n";
       }
       else
       {
	print "HERE\n";
         `$com`;
       }
     }
     else
     {
       if (@files>1)
       {
        `cp $files[0] $outdir/$libname\_1.notR`;
        `cp $files[1] $outdir/$libname\_2.notR`;  
       }
       else
       { 
        `cp $files[0] $outdir/$libname.notR`;
       }
      }
    }
  }
}

sub mapMultipleFiles
{
    my ($files, $jobsubmit, $servicename, $bowtiecmd, $ribosomeInd, $outdir) = @_;
    my $libname=@{$files}[1];
    my @f1=<${$files}[2]>;
    my @f2=();
    if (${$files}[3]!~/^$/)
    {
      @f2=<@{$files}[3]>;
    }
    for (my $i=0; $i<@f1; $i++)
    {
	my $str_file="";
	if (defined $f1[$i])
        {
	  $str_file=$f1[$i]; 
	  if (defined $f2[$i])
	  {
	    $str_file="-1 ".$f1[$i]." -2 ".$f2[$i];
	  }
	  
	  my $com="$bowtiecmd --un-conc $outdir/${libname}_".($i+1).".notR -x $ribosomeInd $str_file --al-conc $outdir/${libname}_".($i+1).".yesR  >  $outdir/${libname}_".($i+1).".mapped";
	  my $job=$jobsubmit." -n ".$servicename."_".$libname."_".($i+1)." -c \"$com\"";
          print $job."\n";   
          `$job`;
	}
    }
}

__END__


=head1 NAME

stepRibo.pl

=head1 SYNOPSIS  

stepRibo.pl -i input <fastq> 
            -o outdir <output directory> 
            -b bowtieCmd <bowtie dir and file> 
            -p params <bowtie params> 
            -r ribosomeInd <ribosome Index file>

stepRibo.pl -help

stepRibo.pl -version

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


stepRibo.pl -i test1.fastq:test2.fastq:ctrl1.fastq:ctrl2.fastq
            -o ~/out
            -b ~/bowtie_dir/bowtie
            -p "-p 8 -n 2 -l 20 -M 1 -a --strata --best"
            -r ~/bowtie_ind/rRNA

=head1 AUTHORS

 Hennady Shulha, PhD 

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


