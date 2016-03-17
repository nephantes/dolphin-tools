#!/usr/bin/env perl

#########################################################################################
#                                      parseBow.pl
#########################################################################################
# 
# Parse bowtie output
#
#########################################################################################
# AUTHORS:
#
# Alper Kucukural, PhD
# 11/03/2015
# 
#########################################################################################


############## LIBRARIES AND PRAGMAS ################

 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 
#################### VARIABLES ######################
 my $name             = "";
 my $file             = "";
 my $paired           = "";
 my $print_version    = "";
 my $help             = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $command=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
    'file=s'         => \$file,
    'name=s'         => \$name,
    'paired=s'       => \$paired,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($name eq "") or ($file eq "")  );	

open(IN, $file);

my $i = 0;
my ($RDS_T, $RDS_P, $RDS_C1, $RDS_C2, $ALGN_T, $a, $b)=(0, 0, 0, 0, 0, 0, 0);
while(my $line=<IN>)
{
  chomp($line);
  $line=~s/^\s+//;
  my @arr=split(/\s/, $line);
  $RDS_T=$arr[0] if ($i=~/^1$/);
  $RDS_P=$arr[0]." ".$arr[1] if ($i == 2);
  
  if ($i == 3)
  {
    $a=$arr[0];
    $RDS_C1=$arr[0]." ".$arr[1]
  }
  if ($i == 4)
  {
    $b=$arr[0];
    $RDS_C2=$arr[0]." ".$arr[1];
  }
  $ALGN_T=($a+$b)." (".$arr[0].")" if (($i == 5 && lc($paired) ne "paired" ) || ($i == 13 && lc($paired) eq "paired" )) ;

  $i++;
}
print "$name\t$RDS_T\t$RDS_P\t$RDS_C1\t$RDS_C2\t$ALGN_T\n";

__END__


=head1 NAME

parseBow.pl

=head1 SYNOPSIS  

parseBow.pl 
            -n name <sample name> 
            -f filename <bowtie std output> 
            -p paired <paired|no> 

parseBow.pl -help

parseBow.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -f filename <output directory>

owtie std output

=head2 -n name <picard running line> 

sample name

=head2  -p <paired>  

paired end library or not <paired|no>

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

This program runs the cufflinks after tophat mappings

=head1 EXAMPLE

parseBow.pl 
            -n name <sample name> 
            -f filename <bowtie std output> 
            -p paired <yes|no> 

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

