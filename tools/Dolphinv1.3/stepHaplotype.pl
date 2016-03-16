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
my $previous         = "";
my $haploCmd         = "";
my $picardCmd        = "";
my $smctfc           = "";
my $smctfe           = "";
my $mbqs             = "";
my $mrpas            = "";
my $mrirps           = "";
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
	'previous=s'     => \$previous,
	'haplocmd=s'     => \$haploCmd,
	'picardCmd=s'    => \$picardCmd,
	'smctfc=s'       => \$smctfc,
    'smctfe=s'       => \$smctfe,
    'mbqs=s'         => \$mbqs,
    'mrpas=s'        => \$mrpas,
    'mrirps=s'       => \$mrirps,
	'jobsubmit=s'    => \$jobsubmit,
	'servicename=s'  => \$servicename,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($genome eq "") or ($outdir eq "") or ($haploCmd eq "") or ($picardCmd eq ""));	

 
################### MAIN PROGRAM ####################
my $inputdir="";
print "$previous\n";

my $tophat_test = `ls $outdir/tophat`;
my $tdf_tophat_test = `ls $outdir/tdf_tophat`;
my $chip_test = `ls $outdir/seqmapping/chip`;
my $input_file_cmd = "";

if ($tophat_test !~/No such file or directory/) {
	$inputdir = "$outdir/tophat";
	$input_file_cmd = "ls $inputdir/*/*.sorted.bam 2>&1";
}elsif ($tdf_tophat_test !~/No such file or directory/) {
	$inputdir = "$outdir/tdf_tophat";
	$input_file_cmd = "ls $inputdir/*/*.sorted.bam 2>&1";
}elsif ($chip_test !~/No such file or directory/){
	$inputdir = "$outdir/seqmapping/chip";
	$input_file_cmd = "ls $inputdir/*.sorted.bam 2>&1";
}else{
	die "Error 256: cannot find bam file directory:";
}

$outdir   = "$outdir/variant_calls";
`mkdir -p $outdir`;
die "Error 15: Cannot create the directory:".$outdir if ($?);

my $com="";
$com=`$input_file_cmd`;
die "Error 64: please check the if you defined the parameters right:" unless ($com !~/No such file or directory/);

print $com;
my @files = split(/[\n\r\s\t,]+/, $com);

foreach my $file (@files)
{
	my $input_file = $file;
	$file=~/.*\/(.*)/;
	my $str_file=$1;
	die "Error 64: please check the file:".$file unless (checkFile($file));
	
	$str_file=~/(.*).*\./;
	my $bname=$1;
	
	############	Start Program Calls
	my $com="";
	##	Make sure bam isn't malformed
	$com.="$picardCmd AddOrReplaceReadGroups I=$input_file O=$outdir/hg_$str_file RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20";
	$com.=" && " if ($com!~/^$/);
	##	Make sure bam is ordered properly
	$com.=" $picardCmd ReorderSam I=$outdir/hg_$str_file O=$outdir/reordered_$str_file R=$genome CREATE_INDEX=TRUE";
	##	Remove hg_file
	$com.=" && rm $outdir/hg_$str_file " if ($com!~/^$/);
	$com.=" && " if ($com!~/^$/);
	##	Run HaplotypeCaller
	$com.="$haploCmd -R $genome -T HaplotypeCaller -I $outdir/reordered_$str_file -o $outdir/$bname.vcf --filter_reads_with_N_cigar --fix_misencoded_quality_scores";
	$com.=" --standard_min_confidence_threshold_for_calling $smctfc --standard_min_confidence_threshold_for_emitting $smctfe --min_base_quality_score $mbqs";
	$com.=" --minReadsPerAlignmentStart $mrpas --maxReadsInRegionPerSample $mrirps";
	##	Remove reordered bam files
	$com.=" && rm $outdir/reordered_$str_file " if ($com!~/^$/);
	$com.=" && rm $outdir/reordered_$bname.bai " if ($com!~/^$/);
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
            -ha haploCmd <haplotypecaller command>
			-pi picardCmd <picard command>

stepHaplotype.pl -help

stepHaplotype.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -o outdir <output directory>

the output files will be "$outdir/variant_calls" 

=head2  -g genome <genome files> 

Genome fasta file. (Full path)

=head2 -ha haploCmd <Haplotypecaller dir and file> 

Full path of haplotypecaller

=head2  -pi picardCmd <picard dir and file>

Full path of picard

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
            -ha haploCmd <haplotypecaller dir and file>
			-pi picardCmd <picard dir and file>

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




