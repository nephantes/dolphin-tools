#!/usr/bin/env perl

#########################################################################################
#                                       stepMergeRSeQC.pl
#########################################################################################
# 
#  This program merges the multiple sample RSeQC output into single file 
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
 use Data::Dumper;
 
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
# runs the RSeQC program

my $indir  = "$outdir/RSeQC_$type";
my $outd  = "$outdir/RSeQC_$type";

`mkdir -p $outd`;
die "Error 15: Cannot create the directory:".$outd if ($?);

my $puboutdir   = "$pubdir/$wkey";
`mkdir -p $puboutdir`;
die "Error 15: Cannot create the directory:".$puboutdir if ($?);

my @files=();
print $type."\n";
my @outtypes = ("RSqQC");
my @order=( "Total Reads", "Total Tags" , "Total Assigned Tags", "CDS_Exons", "5'UTR_Exons", "3'UTR_Exons", "Introns", "TSS_up_1kb", "TSS_up_5kb", "TSS_up_10kb", "TES_down_1kb", "TES_down_5kb", "TES_down_10kb");
my %lines=(
  "Total Reads" => 1,
  "Total Tags" => 1,
  "Total Assigned Tags" => 1,
  "CDS_Exons" => 2,
  "5'UTR_Exons" => 2,
  "3'UTR_Exons" => 2,
  "Introns" => 2,
  "TSS_up_1kb" => 2,
  "TSS_up_5kb" => 2,
  "TSS_up_10kb" => 2,
  "TES_down_1kb" => 2,
  "TES_down_5kb" => 2,
  "TES_down_10kb" => 2
);

my $c=0;
foreach my $outtype (@outtypes)
{

my $ext=".out";
@files = <$indir/$outtype*$ext>;

my @rowheaders=();
my @libs=();
my %vals=();
my %normvals=();
$version="RSeQC.v2.6.2";

foreach my $d (@files){
  my $libname=basename($d, $ext);
  $libname=~s/RSqQC.//g;
  $libname=~s/rsem.out.//g;
  $libname=~s/.genome//g;
  print $libname."\n";
  push(@libs, $libname); 
  getVals($d, $libname, \%vals, \%normvals, \%lines);
}
#print Dumper(%vals);
#print Dumper(%normvals);


my $sizemetrics = keys %vals;
write_results("$outd/$outtype.$type.counts.tsv", \@libs,\%vals, \@order, "region") if ($sizemetrics>0);
write_results("$outd/$outtype.$type.tagskb.tsv", \@libs,\%normvals, \@order, "region") if ($sizemetrics>0);

#Copy count directory to its web accessible area

#my $com="rm -rf $outd/*.$outtype.out && cp -R $outd $puboutdir && ";  
my $com="cp -R $outd $puboutdir ";  

$com.= "&& echo \"$wkey\t$version\tRSeQC_$type\tRSeQC_$type/$outtype.$type.counts.tsv\" >> $puboutdir/reports.tsv "; 
$com.= "&& echo \"$wkey\t$version\tRSeQC_$type\tRSeQC_$type/$outtype.$type.tagskb.tsv\" >> $puboutdir/reports.tsv "; 

print $com."\n"; 
`$com`;
$c++;
}

sub write_results
{
  my ($outfile, $libs, $vals, $order, $name )=@_;
  open(OUT, ">$outfile");
  print OUT "$name\t".join("\t", @{$libs})."\n";

  my $lib=${$libs}[0];
  foreach my $key ( @order )
  {
    if (exists ${$vals}{$lib}{$key}) {
    print OUT $key;
    foreach my $lib (@{$libs})
    {
      print OUT "\t".${$vals}{$lib}{$key};
    } 
    print OUT "\n";
    }
  }
  close(OUT);
}

sub getVals{
  my ($filename, $libname, $vals, $normvals, $lines)=@_;
  if (-e $filename){
     open(IN, $filename);
     while(my $line=<IN>)
     {
       chomp($line);
       my @vals_arr=split(/\s{2,}/,$line);
       if (exists ${$lines}{$vals_arr[0]}) {
         my $idx=${$lines}{$vals_arr[0]};
         ${$vals}{$libname}{$vals_arr[0]}=$vals_arr[$idx] if (exists $vals_arr[$idx]);
         if ($idx==2) {
             ${$normvals}{$libname}{$vals_arr[0]}=$vals_arr[3] if (exists $vals_arr[3]);
         }
       }
     } 
  }
  
}

__END__


=head1 NAME

stepMergeRSeQC.pl

=head1 SYNOPSIS  

stepMergeRSeQC.pl 
            -o outdir <output directory> 
            -t type <Tophat|RNASeq>
            -p pubdir <the path of publicly accessible dir>
            -w wkey <key of the run> 

stepMergeRSeQC.pl -help

stepMergeRSeQC.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -o outdir <output directory>

the output files will be stored "$outdir/RSeQC_$type" 

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

This program merges the multiple sample RSeQC output into single file 

=head1 EXAMPLE

stepMergeRSeQC.pl 
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

