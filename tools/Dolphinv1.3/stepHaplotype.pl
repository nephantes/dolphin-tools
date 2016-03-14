#!/usr/bin/env perl

#########################################################################################
#                                       stepTophat2.pl
#########################################################################################
# 
#  This program uses GATK Haplotype Caller for variant calling 
#
#########################################################################################
# AUTHORS:
#
# Nicholas Merowsky
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
my $genome           = "";
my $paramshaplo     = "";
my $previous         = "";
my $haploCmd         = "";
my $picardCmd        = "";
my $jobsubmit        = "";
my $servicename      = "";
my $help             = "";
my $print_version    = "";
my $version          = "1.0.0";

################### PARAMETER PARSING ####################

my $command=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
	'outdir=s'       => \$outdir,
	'genome=s'       => \$genome,
	'haplocmd=s'     => \$haploCmd,
	'paramshaplo=s'  => \$paramshaplo,
	'previous=s'     => \$previous,
	'jobsubmit=s'    => \$jobsubmit,
	'servicename=s'  => \$servicename,
	'picardCmd=s'    => \$picardCmd,
	'help'           => \$help, 
	'version'        => \$print_version,
) or die("Unrecognized options.\nFor help, run this script with -help option.\n");

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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($genome eq "") or ($outdir eq "") or ($haploCmd eq "") );	

 
################### MAIN PROGRAM ####################
my $inputdir="";
print "$previous\n";
my $tophat_test = `ls $outdir/tdf_tophat`;
if ($tophat_test !~/No such file or directory/) {
	$inputdir = "$outdir/tdf_tophat";
}else{
	my $chip_test = `ls $outdir/macs`;
	if ($chip_test !~/No such file or directory/) {
		$inputdir = "$outdir/macs";
	}
	die "Error 256: cannot find bam file directory:" unless ($chip_test !~/No such file or directory/);
}

$outdir   = "$outdir/variant_calls";
`mkdir -p $outdir`;
die "Error 15: Cannot create the directory:".$outdir if ($?);

$paramshaplo =~s/,/ /g;
$paramshaplo=~s/_/--/g;
my $com="";
$com=`ls $inputdir/*.bam 2>&1`;
die "Error 64: please check the if you defined the parameters right:" unless ($com !~/No such file or directory/);

print $com;
my @files = split(/[\n\r\s\t,]+/, $com);

if (lc($paramshaplo)=~/^no/ || $paramshaplo=~/^$/)
{
	$paramshaplo="";
}
foreach my $file (@files)
{ 
	$file=~/.*\/(.*).bam/;
	my $bname=$1;
	my $str_files=$file;
	die "Error 64: please check the file:".$file unless (checkFile($file));
	
	############	Start Program Calls
	my $com="";
	##	Make sure bam isn't malformed
	$com.="$picardCmd AddOrReplaceReadGroups I=$inputdir/$file O=$outdir/hg_$file RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20";
	$com.=" && " if ($com!~/^$/);
	##	Make sure bam is ordered properly
	$com.=" $picardCmd ReorderSam I=$outdir/hg_$file O=$outdir/reordered_$file R=$genome CREATE_INDEX=TRUE";
	$com.=" && rm $outdir/hg_$file " if ($com!~/^$/);
	$com.=" && " if ($com!~/^$/);
	##	Run HaplotypeCaller
	$com="$haploCmd -R $genome -T HaplotypeCaller -I $outdir/reordered_$file -O $outdir/$file.vcf $paramshaplo --filter_reads_with_N_cigar --fix_misencoded_quality_scores";
	############	End Program Calls
	
	my $job=$jobsubmit." -n ".$servicename."_".$bname." -c \"$com\"";
	print $job."\n";
	`$job`;
	die "Error 25: Cannot run the job:".$job if ($?);
}

sub checkFile
{
	my ($file) = $_[0];
	return 1 if (-e $file);
	return 0;
}

__END__


=head1 NAME

stepHaplotype.pl

=head1 SYNOPSIS  

stepHaplotype.pl 
            -o outdir <output directory> 
            -g genome <genome fasta file> 
            -c cmdhaplo <haplotypecaller Command> 
            -p paramshaplo <haplotypecaller parameters>

stepHaplotype.pl -help

stepHaplotype.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -o outdir <output directory>

the output files will be "$outdir/variant_calls" 

=head2  -g genome <genome files> 

Genome fasta file. (Full path)

=head2 -c CmdHaplo <Haplotypecaller dir and file> 

Fullpath of tophat-calculate-expression file. Ex: /isilon_temp/garber/bin/RSEM/tophat-calculate-expression

=head2  -p haploparams <haplotypecaller additional params>

tophat reference file

=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

This program finds SNP variants using GATK HaplotypeCaller

=head1 EXAMPLE

stepHaplotype.pl 
            -o outdir <output directory> 
            -g genome <genome fasta file> 
            -c haploCmd <haplotypecaller dir and file> 
            -p paramshaplo <haplotypecaller additional params>

=head1 AUTHORS

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




