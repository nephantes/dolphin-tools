#!/usr/bin/env perl

#########################################################################################
#                                       stepMergePicard.pl
#########################################################################################
# 
#  This program merges the multiple sample picard output into single file 
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
 my $type             = "";
 my $mergepicard      = "";
 my $pubdir           = "";
 my $wkey             = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions(
    'outdir=s'        => \$outdir,
    'type=s'          => \$type,
    'mergepicard=s'   => \$mergepicard,
    'pubdir=s'        => \$pubdir,
    'wkey=s'          => \$wkey,
    'help'            => \$help, 
    'version'         => \$print_version,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($type eq "") or ($outdir eq "") );	

################### MAIN PROGRAM ####################
# runs the picard program

my $indir  = "$outdir/picard_$type";
my $outd  = "$outdir/picard_$type";


`mkdir -p $outd`;
die "Error 15: Cannot create the directory:".$outd if ($?);

my $puboutdir   = "$pubdir/$wkey";
`mkdir -p $puboutdir`;
die "Error 15: Cannot create the directory:".$puboutdir if ($?);

my @files=();
print $type."\n";
my @outtypes = ("CollectRnaSeqMetrics", "alignment_summary_metrics", "base_distribution_by_cycle_metrics", "insert_size_metrics", "quality_by_cycle_metrics", "quality_distribution_metrics" );

my $c=0;
foreach my $outtype (@outtypes)
{
my $ext="_multiple.out";
$ext.=".$outtype" if ($outtype ne "CollectRnaSeqMetrics");
print $ext."\n";
@files = <$indir/*$ext>;

my @rowheaders=();
my @libs=();
my %metricvals=();
my %histvals=();
$version="picard_tools_1.131";

my $pdffile="";
foreach my $d (@files){
  my $libname=basename($d, $ext);
  print $libname."\n";
  push(@libs, $libname); 
  getMetricVals($d, $libname, \%metricvals, \%histvals,\@rowheaders);


  if (-e "$outd/".$libname."_multi")
  {
    $pdffile.= "&& $mergepicard $outd/".$libname."_multi/*.pdf $outd/".$libname."_multi_metrics.pdf && rm -rf $outd/".$libname."_multi/";
    $pdffile.= "&& echo \"$wkey\t$version\tpicard_$type\tpicard_$type/".$libname."_multi_metrics.pdf\" >> $puboutdir/reports.tsv "; 
  }
}

my $sizemetrics = keys %metricvals;
write_results("$outd/picard.$outtype.stats.tsv", \@libs,\%metricvals, \@rowheaders, "metric") if ($sizemetrics>0);

my $sizehist = keys %histvals;
write_results("$outd/picard.$outtype.hist.tsv", \@libs,\%histvals, "none", "nt") if ($sizehist>0);



#Copy count directory to its web accessible area

#my $com="rm -rf $outd/*.$outtype.out && cp -R $outd $puboutdir && ";  
my $com="cp -R $outd $puboutdir ";  

$com.= "&& echo \"$wkey\t$version\tpicard_$type\tpicard_$type/picard.$outtype.stats.tsv\" >> $puboutdir/reports.tsv " if ($sizemetrics>0); 
$com.= "&& echo \"$wkey\t$version\tpicard_$type\tpicard_$type/picard.$outtype.hist.tsv\" >> $puboutdir/reports.tsv " if ($sizehist>0); 
$com.= $pdffile;
print $com."\n"; 
`$com`;
$c++;
}

sub write_results
{
  my ($outfile, $libs, $vals, $rowheaders, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\t".join("\t", @{$libs})."\n";
  my $size=0;
  $size=scalar(@{${$vals}{${$libs}[0]}}) if(exists ${$libs}[0] and exists ${$vals}{${$libs}[0]} );

  for (my $i=0; $i<$size;$i++)
  { 
    my $rowname=$i;
    $rowname = ${$rowheaders}[$i] if ($name=~/metric/);
    print OUT $rowname;
    foreach my $lib (@{$libs})
    {
      print OUT "\t".${${$vals}{$lib}}[$i];
    } 
    print OUT "\n";
  }
  close(OUT);
}

sub getMetricVals{
  my ($filename, $libname, $metricvals, $histvals,$rowheaders)=@_;
  if (-e $filename){
     my $nextisheader=0;
     my $nextisvals=0;
     my $nexthist=0;
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       @{$rowheaders}=split(/\t/, $line) if ($nextisheader && !scalar(@{$rowheaders})); 
       if ($nextisvals) {
         @{${$metricvals}{$libname}}=split(/\t/, $line);
         $nextisvals=0;
       }
       if($nexthist){
          my @vals=split(/[\s\t]+/,$line); 
          push(@{${$histvals}{$libname}}, $vals[1]) if (exists $vals[1]);
       }
       $nextisvals=1 if ($nextisheader); $nextisheader=0;
       $nextisheader=1 if ($line=~/METRICS CLASS/);
       $nexthist=1 if ($line=~/normalized_position/);
     } 
  }
  
}

__END__


=head1 NAME

stepMergePicard.pl

=head1 SYNOPSIS  

stepMergePicard.pl 
            -o outdir <output directory> 
            -t type <Tophat|ChipSeq>
            -p pubdir <the path of publicly accessible dir>
            -w wkey <key of the run> 

stepMergePicard.pl -help

stepMergePicard.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -o outdir <output directory>

the output files will be stored "$outdir/picard_$type" 

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

This program merges the multiple sample picard output into single file 

=head1 EXAMPLE

stepMergePicard.pl 
            -o outdir <output directory> 
            -t type <Tophat|ChipSeq>
            -p pubdir <the path of publicly accessible dir>
            -w wkey <key of the run> 

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

