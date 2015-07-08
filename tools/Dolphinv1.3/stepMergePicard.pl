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

@files = <$indir/*.out>;

my @rowheaders=();
my @libs=();
my %metricvals=();
my %histvals=();

foreach my $d (@files){ 
  my $libname=basename($d, ".out");
  push(@libs, $libname); 
  getMetricVals($d, $libname, \%metricvals, \%histvals,\@rowheaders);
}

write_results("$outd/picard.stats.tsv", \@libs,\%metricvals, \@rowheaders, "metric");
write_results("$outd/picard.hist.tsv", \@libs,\%histvals, "none", "nt");

#Copy count directory to its web accessible area

my $com="rm -rf $outd/*.out && cp -R $outd $puboutdir && ";  
$version="picard_tools_1.131";
$com.= "echo \"$wkey\t$version\tpicard.stats\tpicard_$type/picard.stats.tsv\" >> $puboutdir/reports.tsv &&"; 
$com.= "echo \"$wkey\t$version\tpicard.hist\tpicard_$type/picard.hist.tsv\" >> $puboutdir/reports.tsv "; 
print $com."\n"; 
`$com`;

sub write_results
{
  my ($outfile, $libs, $vals, $rowheaders, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\t".join("\t", @{$libs})."\n";
  my $size=0;
  $size=scalar(@{${$vals}{${$libs}[0]}}) if(exists ${$libs}[0]);

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
       @{$rowheaders}=split(/[\s\t]+/, $line) if ($nextisheader && !scalar(@{$rowheaders})); 
       if ($nextisvals) {
         @{${$metricvals}{$libname}}=split(/[\s\t]+/, $line);
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

