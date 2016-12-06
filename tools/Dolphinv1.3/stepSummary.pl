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
 my $pubdir           = "";
 my $wkey             = "";
 my $samtools         = "";
 my $username         = "";
 my $config           = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $command=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
	'outdir=s'       => \$outdir,
    'pubdir=s'       => \$pubdir,
    'wkey=s'         => \$wkey,
	'samtools=s'     => \$samtools,
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

my $reportfile   = "$pubdir/$wkey/reports.tsv";
my $reportsummarydir = "$pubdir/$wkey/summary";
my $outd   = "$outdir/summary";
`mkdir -p $reportsummarydir`;
`mkdir -p $outd`;

my %tsv;
my @headers = ();
my $count_files = `ls $outdir/counts/*.summary.tsv 2>/dev/null`;

push(@headers, 'Sample');
push(@headers, 'Total Reads');

if ( $count_files ne ""){
	my @files = split(/[\n\r\s\t,]+/, $count_files);
	my $filestr="";
	foreach my $file (@files)
	{
		if ($file eq $files[-1]) {
			#append final counts
			parseRNACountFile($file, 1);
		}else{
			#only add specific counts
			parseRNACountFile($file, 0);
		}
	}
}else{
	parseNormalCounts();
}

my $rsem_dir = getDirectory($outdir, 'rsem');
if ($rsem_dir ne "") {
	checkAlignmentType($rsem_dir, "rsem");
}

my $tophat_dir = getDirectory($outdir, 'tophat');
if ($tophat_dir ne "") {
	checkAlignmentType($tophat_dir, "tophat");
}

my $star_dir = getDirectory($outdir, 'star');
if ($star_dir ne "") {
	checkAlignmentType($star_dir, "star");
}

my $hisat2_dir = getDirectory($outdir, 'hisat2');
if ($hisat2_dir ne "") {
	checkAlignmentType($hisat2_dir, "hisat2");
}

my $chip_dir = getDirectory($outdir, 'chip');
if ($chip_dir ne "") {
	checkAlignmentType($chip_dir, "chip");
}else{
	my $chip_seq_dir = getDirectory($outdir, "seqmapping/chip");
	if ($chip_seq_dir ne "") {
		checkAlignmentType($chip_seq_dir, "chip");
	}
}

my $atac_dir = getDirectory($outdir, 'atac');
if ($atac_dir ne "") {
	checkAlignmentType($atac_dir, "atac");
}else{
	my $atac_seq_dir = getDirectory($outdir, "seqmapping/atac");
	if ($atac_seq_dir ne "") {
		checkAlignmentType($atac_seq_dir, "atac");
	}
}

my @keys = keys %tsv;
my $summary = "$outd/summary_data.tsv";
my $header_string = join("\t", @headers);
`echo "$header_string" > $summary`;
`echo "$header_string" > $reportsummarydir/summary.tsv`;
foreach my $key (@keys){
	my $values = join("\t", @{ $tsv{$key} });
	`echo "$values" >> $summary`;
	`echo "$values" >> $reportsummarydir/summary.tsv`;
}

`echo "$wkey\t$version\tsummary\tsummary/summary.tsv" >> $reportfile`;

###############
# SUBROUTINES #
###############

sub parseRNACountFile
{
	my ($file) = $_[0];
	my ($end_file) = $_[1];
	my $contents_full = `cat $file`;
	my @contents_array = split(/[\n\r,]+/, $contents_full);
	my %unmapped_dir = {};
	foreach my $contents_sample (@contents_array)
	{
		my @contents_sample_array = split(/[\t,]+/, $contents_sample);
		if ($contents_sample_array[0] ne 'File') {
			if ($tsv{$contents_sample_array[0]} eq undef) {
				$tsv{$contents_sample_array[0]} = [$contents_sample_array[0], $contents_sample_array[1]];
			}
			if ($tsv{$contents_sample_array[0]}[1] != undef) {
				if ($tsv{$contents_sample_array[0]}[1] < $contents_sample_array[1]) {
					$tsv{$contents_sample_array[0]}[1]= $contents_sample_array[1];
				}
			}
			if ($tsv{$contents_sample_array[0]}[1] != undef) {
				if ($tsv{$contents_sample_array[0]}[1] < $contents_sample_array[1]) {
					$tsv{$contents_sample_array[0]}[1]= $contents_sample_array[1];
				}
			}
			
			my @totalreads = split(/[\s\t,]+/, $contents_sample_array[2]);
			if ($unmapped_dir{$contents_sample_array[0]} == undef || $unmapped_dir{$contents_sample_array[0]} >= $totalreads[0]) {
				$unmapped_dir{$contents_sample_array[0]} = $totalreads[0];
			}
			
			my @reads1 = split(/[\s,]+/, $contents_sample_array[5]);
			push($tsv{$contents_sample_array[0]}, $reads1[0]);
			if ($end_file) {
				push($tsv{$contents_sample_array[0]}, $unmapped_dir{$contents_sample_array[0]});
			}
		}
	}
	my @split_file = split(/[\/]+/, $file);
	my $name = $split_file[-1];
	$name =~/(.*)\.summary\.tsv/;
	my $split_name = $1;
	push(@headers, $split_name);
	if ($end_file) {
		push(@headers, 'Total align');
	}
}

sub parseNormalCounts
{
	chomp(my $contents = `ls $outdir/input`);
	my @files = split(/[\n]+/, $contents);
	foreach my $file (@files){
		$file=~/(.*)\..*/;
		my $split_name = $1;
		$split_name=~s/\.1$//g;
		$split_name=~s/\.2$//g;
		print $split_name;
		chomp(my $rawcount = `wc -l $outdir/input/$file`);
		my @counts = split(/[ ]+/, $rawcount);
		my $finalcount = $counts[0]/4;
		$tsv{$split_name} = [$split_name, $finalcount];
	}
}

sub checkAlignmentType
{
	my $directories = $_[0];
	my $type = $_[1];
	my $deduptype = $type;
	$deduptype = $type."_ref.transcripts" if ($type eq "rsem");
	my @dirs = split(/[\n]+/, $directories);
	if(grep( /^$outdir\/dedupmerge$deduptype$/, @dirs )) {
		dedupReadsAligned("$outdir/dedupmerge$deduptype", $type, "merge");
	}elsif(grep( /^$outdir\/dedup$deduptype$/, @dirs )){
		dedupReadsAligned("$outdir/dedup$deduptype", $type, "norm");
	}elsif(grep( /^$outdir\/merge$type$/, @dirs )){
		searchAligned("$outdir/merge$type", $type, "*.bam", "merge");
	}elsif(grep( /^$outdir\/$type$/, @dirs )){
		if ($type eq "tophat"){
			alteredAligned("$outdir/$type", $type, "*/*.sorted.bam", "norm");
		}elsif ($type eq "rsem"){
			my $genome_check = `ls $outdir/$type/*/*genome.bam 2>&1`;
			if ($genome_check !~ /No such file or directory/) {
				alteredAligned("$outdir/$type", $type, "*/*genome.bam", "norm");
			}else{
				alteredAligned("$outdir/$type", $type, "*/*transcript.bam", "norm");
			}
		}else{
			alteredAligned("$outdir/$type", $type, "*/*.sorted.bam", "norm");
		}
	}elsif($type eq "chip"){
		searchAligned("$outdir/seqmapping/chip", $type, "*.bam", "norm");
	}elsif($type eq "atac"){
		searchAligned("$outdir/seqmapping/atac", $type, "*.sorted.b*", "norm");
	}
}

sub dedupReadsAligned
{
	my ($directory) = $_[0];
	my ($type) = $_[1];
	my ($merge) = $_[2];
	chomp(my $contents = `ls $directory/*PCR_duplicates`);
	my @files = split(/[\n]+/, $contents);
	push(@headers, "Duplicated Reads $type");
	push(@headers, "Multimapped Reads Aligned $type");
	push(@headers, "Unique Reads Aligned $type");
	foreach my $file (@files){
		my $multimapped;
		$file=~/.*\/(.*)/;
		my $file_solo = $1;
		my @nodup = split(/_PCR_duplicates/, $file_solo);
		my $name = $nodup[0];
		if ($type eq 'rsem') {
			print "awk 'NR == 2 {print \$3}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt \n";
			chomp($multimapped = `awk 'NR == 2 {print \$3}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt`)
		}elsif($type eq "tophat"){
			print "cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep -A 1 'Aligned pairs:' | awk 'NR % 3 == 2 {sum+=\$3} END {print sum}' \n";
			chomp($multimapped = `cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep -A 1 'Aligned pairs:' | awk 'NR % 3 == 2 {sum+=\$3} END {print sum}'`)
		}elsif($type eq "chip"){
			print "cat $outdir/seqmapping/chip/$name*.sum | awk '{sum+=\$7} END {print sum}' \n";
			chomp($multimapped = `cat $outdir/seqmapping/chip/$name*.sum | awk '{sum+=\$7} END {print sum}'`);
		}elsif($type eq "atac"){
			print "cat $outdir/seqmapping/atac/$name*.sum | awk '{sum+=\$7} END {print sum}' \n";
			chomp($multimapped = `cat $outdir/seqmapping/atac/$name*.sum | awk '{sum+=\$7} END {print sum}'`);
		}else{
			print "$samtools view -f 256 $directory/$name*.bam | awk '{print \$1}' | sort -u | wc -l \n";
			chomp($multimapped = `$samtools view -f 256 $directory/$name*.bam | awk '{print \$1}' | sort -u | wc -l`);
		}
		print "cat $file | grep -A 1 \"LIB\" | grep -v \"LIB\"";
		chomp(my $aligned = `cat $file | grep -A 1 \"LIB\" | grep -v \"LIB\"`);
		my @values = split("\t", $aligned);
		my $dedup = int($values[5]) + int($values[4]);
		my $total = int($values[2]) + int($values[1]) - $dedup - int($multimapped);
		push($tsv{$name}, $dedup.'');
		push($tsv{$name}, $multimapped.'');
		push($tsv{$name}, $total.'');
	}
}

sub searchAligned
{
	my ($directory) = $_[0];
	my ($type) = $_[1];
	my ($filetype) = $_[2];
	my ($merge) = $_[3];
	chomp(my $contents = `ls $directory/$filetype`);
	my @files = split(/[\n]+/, $contents);
	push(@headers, "Multimapped Reads Aligned $type");
	push(@headers, "Unique Reads Aligned $type");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		my @split_name = split(/[\/]+/, $file);
		my @namelist = split(/\.bam/, $split_name[-1]);
		my $name = $namelist[0];
		my @sorted = split(/\.sorted/,$namelist[0]);
		$name = $sorted[0];
		if ($type eq 'rsem') {
			print "awk 'NR == 1 {print \$2}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt \n";
			chomp($aligned = `awk 'NR == 1 {print \$2}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt`);
			print "awk 'NR == 2 {print \$3}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt \n";
			chomp($multimapped = `awk 'NR == 2 {print \$3}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt`);
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, (int($aligned) - int($multimapped))."");
		}elsif($type eq "tophat"){
			my $tophat_parse = "cat $outdir/tophat/pipe.tophat.$name/align_summary.txt";
			if ($merge eq "merge") {
				my $merged_command = "";
				chomp(my $unmerged = `ls -d $outdir/tophat/pipe.tophat.$name*`);
				my @unmerged_dirs = split(/[\n]+/, $unmerged);
				foreach my $unmerge_dir (@unmerged_dirs){
					chomp($unmerge_dir);
					if ($unmerge_dir =~ /$name\_[\d][\d]$/) {
						$merged_command .= "cat" if ($merged_command eq "");
						$merged_command .= " $unmerge_dir/align_summary.txt";
					}
				}
				$tophat_parse = $merged_command;
			}
			print "$tophat_parse | grep 'Aligned pairs:' | awk '{sum+=\$3} END {print sum}' \n";
			chomp($aligned = `$tophat_parse | grep 'Aligned pairs:' | awk '{sum+=\$3} END {print sum}'`);
			if ($aligned eq "") {
				print "$tophat_parse | grep 'Mapped' | awk '{sum+=\$3} END {print sum}' \n";
				chomp($aligned = `$tophat_parse | grep 'Mapped' | awk '{sum+=\$3} END {print sum}'`);
				print "$tophat_parse | grep 'multiple alignments' | awk '{sum+=\$3} END {print sum}' \n";
				chomp($multimapped = `$tophat_parse | grep 'multiple alignments' | awk '{sum+=\$3} END {print sum}'`);
			}else{
				print "$tophat_parse | grep -A 1 'Aligned pairs:' | awk 'NR % 3 == 2 {sum+=\$3} END {print sum}' \n";
				chomp($multimapped = `$tophat_parse | grep -A 1 'Aligned pairs:' | awk 'NR % 3 == 2 {sum+=\$3} END {print sum}'`);
			}
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, (int($aligned) - int($multimapped))."");
		}elsif($type eq "chip"){
			my $chip_parse = "cat $outdir/seqmapping/chip/$name.sum";
			if ($merge eq "merge") {
				my $merged_command = "";
				chomp(my $unmerged = `ls -d $outdir/seqmapping/chip/$name*.sum`);
				my @unmerged_dirs = split(/[\n]+/, $unmerged);
				foreach my $unmerge_dir (@unmerged_dirs){
					chomp($unmerge_dir);
					if ($unmerge_dir =~ /$name\_[\d][\d].sum$/) {
						$merged_command .= "cat" if ($merged_command eq "");
						$merged_command .= " $unmerge_dir";
					}
				}
				$chip_parse = $merged_command;
			}
			print "$chip_parse | awk '{sum+=\$5} END {print sum}' \n";
			chomp($aligned = `$chip_parse | awk '{sum+=\$5} END {print sum}'`);
			print "$chip_parse | awk '{sum+=\$7} END {print sum}' \n";
			chomp($multimapped = `$chip_parse | awk '{sum+=\$7} END {print sum}'`);
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, $aligned);
		}elsif($type eq "atac"){
			my $atac_parse = "cat $outdir/seqmapping/atac/$name.sum";
			if ($merge eq "merge") {
				my $merged_command = "";
				chomp(my $unmerged = `ls -d $outdir/seqmapping/atac/$name*.sum`);
				my @unmerged_dirs = split(/[\n]+/, $unmerged);
				foreach my $unmerge_dir (@unmerged_dirs){
					chomp($unmerge_dir);
					if ($unmerge_dir =~ /$name\_[\d][\d].sum$/) {
						$merged_command .= "cat" if ($merged_command eq "");
						$merged_command .= " $unmerge_dir";
					}
				}
				$atac_parse = $merged_command;
			}
			print "$atac_parse | awk '{sum+=\$5} END {print sum}' \n";
			chomp($aligned = `$atac_parse | awk '{sum+=\$5} END {print sum}'`);
			print "$atac_parse | awk '{sum+=\$7} END {print sum}' \n";
			chomp($multimapped = `$atac_parse | awk '{sum+=\$7} END {print sum}'`);
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, $aligned);
		}else{
			print "$samtools flagstat $file \n";
			chomp($aligned = `$samtools flagstat $file`);
			my @aligned_split = split(/[\n]+/, $aligned);
			my @paired = split(/[\s]+/, $aligned_split[9]);
			my @singleton = split(/[\s]+/, $aligned_split[10]);
			print "$samtools view -f 256 $directory/$name*.bam | awk '{print \$1}' | sort -u | wc -l \n";
			chomp($multimapped = `$samtools view -f 256 $directory/$name*.bam | awk '{print \$1}' | sort -u | wc -l`);
			if (int($paired[0])/2 == 0) {
				print "$samtools view -F 4 $file | awk '{print \$1}' | sort -u | wc -l \n";
				chomp(my $aligned = `$samtools view -F 4 $file | awk '{print \$1}' | sort -u | wc -l`);
				push($tsv{$name}, $multimapped);
				push($tsv{$name}, (int($aligned) - int($multimapped))."");
			}else{
				push($tsv{$name}, $multimapped);
				push($tsv{$name}, ((int($paired[0])/2 - int($singleton[0])) - int($multimapped))."");
			}
		}
	}
}

sub alteredAligned
{
	my ($directory) = $_[0];
	my ($type) = $_[1];
	my ($filetype) = $_[2];
	my ($merge) = $_[3];
	chomp(my $contents = `ls $directory/$filetype`);
	my @files = split(/[\n]+/, $contents);
	push(@headers, "Multimapped Reads Aligned $type");
	push(@headers, "Unique Reads Aligned $type");
	foreach my $file (@files){
		my $multimapped;
		my $aligned;
		$file=~/.*\/(.*)\..*/;
		my $name = $1;
		if ($type eq 'rsem') {
			$name=~s/rsem\.out\.//g;
			$name=~s/\.genome$//g;
			print "awk 'NR == 1 {print \$2}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt \n";
			chomp($aligned = `awk 'NR == 1 {print \$2}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt`);
			print "awk 'NR == 2 {print \$3}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt \n";
			chomp($multimapped = `awk 'NR == 2 {print \$3}' $outdir/rsem/pipe.rsem.$name/rsem.out.$name.stat/rsem.out.$name.cnt`);
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, (int($aligned) - int($multimapped))."");
		}elsif($type eq "tophat"){
			$name=~s/\.sorted$//g;
			print "cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep 'Aligned pairs:' | awk '{sum+=\$3} END {print sum}' \n";
			chomp($aligned = `cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep 'Aligned pairs:' | awk '{sum+=\$3} END {print sum}'`);
			if ($aligned eq "") {
				print "cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep 'Mapped' | awk '{sum+=\$3} END {print sum}' \n";
				chomp($aligned = `cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep 'Mapped' | awk '{sum+=\$3} END {print sum}'`);
				print "cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep 'multiple alignments' | awk '{sum+=\$3} END {print sum}' \n";
				chomp($multimapped = `cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep 'multiple alignments' | awk '{sum+=\$3} END {print sum}'`);
			}else{
				print "cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep -A 1 'Aligned pairs:' | awk 'NR % 3 == 2 {sum+=\$3} END {print sum}' \n";
				chomp($multimapped = `cat $outdir/tophat/pipe.tophat.$name*/align_summary.txt | grep -A 1 'Aligned pairs:' | awk 'NR % 3 == 2 {sum+=\$3} END {print sum}'`);
			}
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, (int($aligned) - int($multimapped))."");
		}elsif($type eq "chip"){
			print "cat $outdir/seqmapping/chip/$name*.sum | awk '{sum+=\$5} END {print sum}' \n";
			chomp($aligned = `cat $outdir/seqmapping/chip/$name*.sum | awk '{sum+=\$5} END {print sum}'`);
			print "cat $outdir/seqmapping/chip/$name*.sum | awk '{sum+=\$7} END {print sum}' \n";
			chomp($multimapped = `cat $outdir/seqmapping/chip/$name*.sum | awk '{sum+=\$7} END {print sum}'`);
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, $aligned);
		}elsif($type eq "atac"){
			print "cat $outdir/seqmapping/atac/$name*.sum | awk '{sum+=\$5} END {print sum}' \n";
			chomp($aligned = `cat $outdir/seqmapping/atac/$name*.sum | awk '{sum+=\$5} END {print sum}'`);
			print "cat $outdir/seqmapping/atac/$name*.sum | awk '{sum+=\$7} END {print sum}' \n";
			chomp($multimapped = `cat $outdir/seqmapping/atac/$name*.sum | awk '{sum+=\$7} END {print sum}'`);
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, $aligned);
		}elsif($type eq "star"){
			$name=~s/\.sorted$//g;
			print "cat $outdir/star/pipe.star.$name*/$name*Log.final.out | grep 'Uniquely mapped reads number' | awk '{sum+=\$6} END {print sum}' \n";
			chomp($aligned = `cat $outdir/star/pipe.star.$name*/$name*Log.final.out | grep 'Uniquely mapped reads number' | awk '{sum+=\$6} END {print sum}'`);
			print "cat $outdir/star/pipe.star.$name*/$name*Log.final.out | grep 'Number of reads mapped to' | awk 'NR % 2 == 1 {sum+=\$9}; NR % 2 == 0 {sum+=\$10} END {print sum}' \n";
			chomp($multimapped = `cat $outdir/star/pipe.star.$name*/$name*Log.final.out | grep 'Number of reads mapped to' | awk 'NR % 2 == 1 {sum+=\$9}; NR % 2 == 0 {sum+=\$10} END {print sum}'`);
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, (int($aligned) - int($multimapped))."");
		}elsif($type eq "hisat2"){
			$name=~s/\.sorted$//g;
			print "cat $outdir/hisat2/pipe.hisat2.$name*/align_summary.txt | grep 'aligned concordantly exactly 1 time' | awk '{sum+=\$1} END {print sum}' \n";
			chomp($aligned = `cat $outdir/hisat2/pipe.hisat2.$name*/align_summary.txt | grep 'aligned concordantly exactly 1 time' | awk '{sum+=\$1} END {print sum}'`);
			print "cat $outdir/hisat2/pipe.hisat2.$name*/align_summary.txt | grep 'aligned concordantly >1 times' | awk '{sum+=\$1} END {print sum}' \n";
			chomp($multimapped = `cat $outdir/hisat2/pipe.hisat2.$name*/align_summary.txt | grep 'aligned concordantly >1 times' | awk '{sum+=\$1} END {print sum}'`);
			push($tsv{$name}, $multimapped);
			push($tsv{$name}, (int($aligned) - int($multimapped))."");
		}else{
			print "$samtools flagstat $file \n";
			chomp($aligned = `$samtools flagstat $file`);
			my @aligned_split = split(/[\n]+/, $aligned);
			my @paired = split(/[\s]+/, $aligned_split[9]);
			my @singleton = split(/[\s]+/, $aligned_split[10]);
			print "$samtools view -f 256 $directory/$name*.bam | awk '{print \$1}' | sort -u | wc -l \n";
			chomp($multimapped = `$samtools view -f 256 $directory/$name*.bam | awk '{print \$1}' | sort -u | wc -l`);
			if (int($paired[0])/2 == 0) {
				print "$samtools view -F 4 $file | awk '{print \$1}' | sort -u | wc -l \n";
				chomp(my $aligned = `$samtools view -F 4 $file | awk '{print \$1}' | sort -u | wc -l`);
				push($tsv{$name}, $multimapped);
				push($tsv{$name}, (int($aligned) - int($multimapped))."");
			}else{
				push($tsv{$name}, $multimapped);
				push($tsv{$name}, ((int($paired[0])/2 - int($singleton[0])) - int($multimapped))."");
			}
		}
	}
}

sub getDirectory
{
	my ($outdir) = $_[0];
	my ($type) = $_[1];
	chomp(my $directories = `ls -d $outdir/*$type* 2>/dev/null`);
	return $directories;
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
