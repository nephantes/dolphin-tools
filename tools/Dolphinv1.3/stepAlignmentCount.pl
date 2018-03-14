#!/usr/bin/env perl

#########################################################################################
#                                       stepSummary.pl
#########################################################################################
# 
#  This program creates a summary file for table generation.
#
#
#########################################################################################
# AUTHORS:
#
# Alper Kucukural, PhD
# Nicholas Merowsky
# Jul 4, 2014
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
 my $samtools         = "";
 my $type             = "";
 my $pubdir           = "";
 my $wkey             = "";
 my $username         = "";
 my $config           = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $command=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
	'outdir=s'       => \$outdir,
	'samtools=s'     => \$samtools,
	'type=s'         => \$type,
    'pubdir=s'       => \$pubdir,
    'wkey=s'         => \$wkey,
    'config=s'       => \$config,
    'username=s'     => \$username,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($outdir eq "") );	

################### MAIN PROGRAM ####################
# Obtain summary information and create file

my $reportfile = "$pubdir/$wkey/reports.tsv";
my $inputdir = "$outdir/$type";
$inputdir = "$outdir/seqmapping/chip" if ($type eq "chip");
$inputdir = "$outdir/seqmapping/atac" if ($type eq "atac");
if ($type eq "tophat" || $type eq "rsem" || $type eq "star" || $type eq "hisat2") {
	my $com=`ls -d $outdir/$type/pipe* 2>&1`;
	my @dirs = split(/[\n\r\s\t,]+/, $com);
	my $bamfile = "accepted_hits.bam";
	$bamfile = "*genome.bam" if ($type eq "rsem");
	$bamfile = "*sorted.bam" if ($type eq "tophat" || $type eq "star" || $type eq "hisat2");
	if ($com !~ /No such file or directory/) {
		foreach my $dir (@dirs)
		{
			if ($bamfile eq "*genome.bam") {
				my $genome_check = `ls $dir/$bamfile 2>&1`;
				if ($genome_check=~/No such file or directory/) {
					$bamfile = "*transcript.bam";
				}
			}
			my @samplename = split(/pipe\.$type\./, $dir);
			my $bname=$samplename[-1];
			$com ="$samtools flagstat $dir/$bamfile > $inputdir/".$bname.".flagstat.txt && ";
			$com.="mkdir -p $pubdir/$wkey/$type && cp $inputdir/".$bname.".flagstat.txt $pubdir/$wkey/$type && ";
			$com.="echo \"$wkey\t$version\tsummary\t$type/$bname.flagstat.txt\" >> $reportfile ";
			my $retval=`$com`;
		}
	}

}else{
	my $com=`ls $inputdir/*.bam 2>&1`;
	my @files = split(/[\n\r\s\t,]+/, $com);
	
	if ($com !~ /No such file or directory/) {
		foreach my $file (@files)
		{
			$file=~/.*\/(.*).bam/;
			my $bname=$1;
			my $bname2 = $bname;
			$bname2=~s/\.sorted//;
			$com ="$samtools flagstat $inputdir/$bname.bam > $inputdir/".$bname2.".flagstat.txt && ";
			$com.="mkdir -p $pubdir/$wkey/$type && cp $inputdir/".$bname2.".flagstat.txt $pubdir/$wkey/$type && ";
			$com.="echo \"$wkey\t$version\tsummary\t$type/$bname2.flagstat.txt\" >> $reportfile ";
			my $retval = `$com`;
		}
	}
}


__END__


=head1 NAME

stepClean.pl

=head1 SYNOPSIS  

stepSummary.pl -o outdir <output directory> 

stepSummary.pl -help

stepSummary.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -o outdir <output directory>

the output files will be "$outdir/after_ribosome" 

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program cleans intermediate files

=head1 EXAMPLE

stepClean.pl 
            -o ~/out

=head1 AUTHORS

 Alper Kucukural, PhD
 Nicholas Merowsky
 
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
