#!/usr/bin/perl

#########################################################################################
#                                       calcCounts.pl
#########################################################################################
# 
#  This program removes adapter sequence. 
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
 my $bamfiles         = "";
 my $bedfile          = "";
 my $outdir           = "";
 my $cmd              = ""; 
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
        'bam=s'          => \$bamfiles,
	'outdir=s'       => \$outdir,
        'bed=s'          => \$bedfile,
        'cmd=s'          => \$cmd,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($bamfiles eq "") or ($outdir eq "") );	

################### MAIN PROGRAM ####################
#    maps the reads to the ribosome and put the files under $outdir/after_ribosome directory


my @bams=split(/[\s,]+/,$bamfiles);
foreach my $bam (@bams)
{
   my $bname=basename($bam, ".bam");
   my $com="$cmd view $bam|awk '{print $3}'|sort|uniq -c|sort -k2,2b>$bname.counts";
   
}


my ( %row_vals, %data );
foreach my $file ( sort glob("*.counts") ) {

    open my $fh, "<", $file or die $!;

    while ( my $line = <$fh> ) {
        chomp $line;
	$line=~s/^[\s\t]+//;
        my ( $values, $row_val ) = split /[\s\t,]+/, $line;
        $row_vals{$row_val} = 1;
        $data{$file}{$row_val} = \$values;
    }
    close $fh;
}
my $i=1;
foreach my $row_val ( sort keys %row_vals ) {
   if($i==1)
   {
      print "id";
      foreach my $file ( sort keys %data ) {
         print "\t".$file; 
      }
   }
   print "\n";
   print $row_val;
   foreach my $file ( sort keys %data ) {
      my %vals=%{$data{$file}};
      if (exists $vals{$row_val}){
         print "\t".$vals{$row_val};
      }
      else{
         print "\t0";
      }
   }
   print "\n";
   $i++;
}


__END__


=head1 NAME

calcCounts.pl

=head1 SYNOPSIS  

calcCounts.pl -i input <fastq> 
            -o outdir <output directory> 
            -b bowtieCmd <bowtie dir and file> 
            -p params <bowtie params> 
            -r ribosomeInd <ribosome Index file>

calcCounts.pl -help

calcCounts.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -b  bam file <bam file> 

multiple vam file can be used

=head2 -o outdir <output directory>

the output files will be "~/out" 

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program quantify frequeincies per region

=head1 EXAMPLE


calcCounts.pl --bam test1.bam,test2.bam
            -o ~/out
            --bed tRNA.bed

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




