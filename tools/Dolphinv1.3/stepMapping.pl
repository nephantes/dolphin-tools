#!/usr/bin/env perl

#########################################################################################
#                                       stepMapping.pl
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
 my $cmd              = "";
 my $spaired          = "";
 my $jobsubmit        = "";
 my $awkdir           = "";
 my $bowtiePar        = "";
 my $advparams        = "";
 my $bowtiecmd        = ""; 
 my $samtoolscmd      = "";
 my $servicename      = "";
 my $param            = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $command=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
    'input=s'        => \$input,
    'outdir=s'       => \$outdir,
    'cmd=s'          => \$bowtiecmd,
    'msamtoolscmd=s' => \$samtoolscmd,
    'awkdir=s'       => \$awkdir,
    'dspaired=s'     => \$spaired,
    'bowtiePar=s'    => \$bowtiePar,
    'servicename=s'  => \$servicename,
    'jobsubmit=s'    => \$jobsubmit,
    'radvparam=s'    => \$advparams,
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
#    maps the reads to the index files and put the files under $outdir/recmapping/$indexname directory

my ($indexfile, $indexname, $indexpar, $description, $filterout, $previous)=split(/,/, $bowtiePar);

$indexpar=~s/__tc__/,/g;
$indexpar=~s/_/ /g;

#if ($advparams ne "NONE")
#{
#  $indexpar=$advparams;
#  $indexpar=~s/,/ /g; 
#  $indexpar=~s/_/-/g; 
#}

print "indexpar=$indexpar\n\n";

my $inputdir="";
print "$previous\n";

if ($previous=~/NONE/)
{
  $inputdir = "$outdir/input";
}
else
{
  $inputdir = "$outdir/seqmapping/".lc($previous);
}
$outdir   = "$outdir/seqmapping/".lc($indexname);


print "I:$inputdir\n";
print "O:$outdir\n";

`mkdir -p $outdir`;

die "Error 15: Cannot create the directory:".$outdir  if ($?);

my $fasta=createIndex($indexfile, $bowtiecmd, $inputdir, $outdir);


my $com="";

if (lc($spaired) =~ /^no/)
{
 $com=`ls $inputdir/*.fastq 2>&1`;
}
else
{
 $com=`ls $inputdir/*.1.fastq 2>&1`;
}
#exit 0 if ($inputdir=~/DC/);
die "Error 64: please check the if you defined the parameters right:" unless ($com !~/No such file or directory/);

print $com;
my @files = split(/[\n\r\s\t,]+/, $com);
$com="";
foreach my $file (@files)
{
  $file=~/.*\/(.*).fastq/;
  my $bname=$1;
  if (lc($spaired) =~ /^no/)
  {
       die "Error 64: please check the file:".$file unless (checkFile($file)); 
       $com="$bowtiecmd $indexpar --no-unal --un $outdir/$bname.fastq -x $indexfile $file --al $outdir/$bname.fastq.mapped -S $outdir/$bname.sam > $outdir/$bname.bow 2>&1 && ";
  }
  else
  {
       $file=~/.*\/(.*).1.fastq/;
       $bname=$1;
       my $str_file="-1 $inputdir/$bname.1.fastq -2 $inputdir/$bname.2.fastq";
       die "Error 64: please check the file: $inputdir/$bname.1.fastq" unless (checkFile("$inputdir/$bname.1.fastq")); 
       die "Error 64: please check the file: $inputdir/$bname.1.fastq" unless (checkFile("$inputdir/$bname.2.fastq")); 
       $com="$bowtiecmd $indexpar --no-unal --un-conc $outdir/$bname.fastq -x $indexfile $str_file --al-conc $outdir/$bname.fastq.mapped -S $outdir/$bname.sam > $outdir/$bname.bow 2>&1 && ";
  }
 $com.="grep -v Warning $outdir/$bname.bow > $outdir/$bname.tmp && ";
 $com.="mv $outdir/$bname.tmp  $outdir/$bname.bow && ";
 $com.="$awkdir/parseBow.pl -n $bname -f $outdir/$bname.bow -p $spaired> $outdir/$bname.sum && ";
 $com.=": \\\$( $samtoolscmd view -bT $fasta $outdir/$bname.sam > $outdir/$bname.tmp1.bam ) && ";
 my $lasttmp="tmp1";
 if (lc($spaired)!~/^no/)
 {
   $com.=": \\\$( $samtoolscmd sort -n $outdir/$bname.tmp1.bam $outdir/$bname.tmp2 ) && "; 
   $com.=": \\\$( $samtoolscmd view  -bf 0x02 $outdir/$bname.tmp2.bam > $outdir/$bname.tmp3.bam ) && "; 
   $lasttmp="tmp3";
 }
 $com.=": \\\$($samtoolscmd sort $outdir/$bname.$lasttmp.bam $outdir/$bname.sorted ) && ";
 $com.=": \\\$($samtoolscmd index $outdir/$bname.sorted.bam ) && ";
 $com.="rm -rf $outdir/$bname.sam && ";
 $com.="rm -rf $outdir/$bname.tmp*.bam && ";
 $com.="rm -rf $outdir/$bname*.mapped";

 my $job=$jobsubmit." -n ".$servicename."_".$bname." -c \"$com\"";
 print "\n\n".$job."\n\n";   
 `$job`;
 die "Error 25: Cannot run the job:".$job if ($?);
}

sub createIndex
{
  my ($ind, $cmd, $inputdir, $outdir) = @_;
     my $fasta=$ind;
     if(checkFile("$ind.fasta"))
     {
       $fasta.=".fasta";
     }
     elsif(checkFile("$ind.fa"))
     {
       $fasta.=".fa";
     }
     else
     {
      $com="ln -s $inputdir/*.fastq $outdir/.";
      `$com`;
      exit(0);
     }
     if(!checkFile("$ind.4.bt2"))
     {
        $com=$cmd."-build $fasta $ind";
        `$com`;
     }
     return $fasta;
}


sub checkFile
{
 my ($file) = $_[0];
 return 1 if (-e $file);
 return 0;
}

__END__


=head1 NAME

stepMapping.pl

=head1 SYNOPSIS  

stepMapping.pl -i input <fastq> 
            -o outdir <output directory> 
            -b bowtieCmd <bowtie dir and file> 
            -p params <bowtie params> 
            -r ribosomeInd <ribosome Index file>

stepMapping.pl -help

stepMapping.pl -version

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


stepMapping.pl -i test1.fastq:test2.fastq:ctrl1.fastq:ctrl2.fastq
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

 
