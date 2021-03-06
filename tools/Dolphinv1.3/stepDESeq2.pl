#!/usr/bin/env perl

#########################################################################################
#                                       stepDESeq2.pl
#########################################################################################
# 
#  This program runs DESeq2 in R.
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
 my $cols             = "";
 my $conds            = "";
 my $num              = "";
 my $heatmap          = "";
 my $fitType          = "";
 my $foldChange       = "";
 my $dataset          = "";
 my $padj             = "";
 my $rscriptCMD       = "";
 my $pubdir           = "";
 my $wkey             = "";
 my $outdir           = "";
 my $jobsubmit        = "";
 my $servicename      = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
    'cols=s'         => \$cols,
    'dconds=s'       => \$conds,
    'dataset=s'      => \$dataset,
    'num=s'          => \$num,
    'eheatmap=s'     => \$heatmap,
    'tfitType=s'     => \$fitType,
    'padj=s'         => \$padj,
    'foldChange=s'   => \$foldChange,
    'rscriptCMD=s'   => \$rscriptCMD,
    'outdir=s'       => \$outdir,
    'pubdir=s'       => \$pubdir,
    'wkey=s'         => \$wkey,
    'servicename=s'  => \$servicename,
    'jobsubmit=s'    => \$jobsubmit,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($cols eq "") or ($conds eq "") or ($outdir eq "") );	

################### MAIN PROGRAM ####################
#    maps the reads to the ribosome and put the files under $outdir/after_ribosome directory

my $inputdir = "$outdir/rsem";
my $input_file_suffix = "_expression_expected_count.tsv";
if (lc($dataset)!~/rsem/)
{
   $inputdir="$outdir/counts";
   $input_file_suffix = ".counts.tsv";
}

$outdir   = "$outdir/DESeq2".$dataset.$num;
`mkdir -p $outdir`;
$cols=~s/[\s\t]+//g;
$conds=~s/[\s\t]+//g;
$cols = checkCols($cols);
$cols=~s/,/\",\"/g;
$cols="c(\"$cols\")";
$cols=~s/-/\./g;
$conds=~s/,/\",\"/g;
$conds="c(\"$conds\")";

my $puboutdir   = "$pubdir/$wkey";
`mkdir -p $puboutdir`;
if (lc($dataset) =~/rsem/ )
{
   makePlot( "genes", $input_file_suffix, $inputdir, $outdir, $cols, $conds, $fitType, $heatmap, $padj, $foldChange, $puboutdir, "DESeq2".$dataset.$num, $wkey);
   makePlot( "isoforms", $input_file_suffix, $inputdir, $outdir, $cols, $conds, $fitType, $heatmap, $padj, $foldChange, $puboutdir, "DESeq2".$dataset.$num, $wkey);
}
else
{
   makePlot( $dataset, $input_file_suffix, $inputdir, $outdir, $cols, $conds, $fitType, $heatmap, $padj, $foldChange, $puboutdir, "DESeq2".$dataset.$num, $wkey);
}

`cp -R $outdir $puboutdir/.`;

sub makePlot
{
my ($type,$input_file_suffix,$inputdir, $outdir, $cols, $conds, $fitType, $heatmap, $padj_cutoff, $foldChange_cutoff, $puboutdir, $deseqdir, $wkey)=@_;
my $output = "$outdir/rscript_$type.R";
my $table = "$outdir/deseq2_$type.tsv";
my $pdfname = "$outdir/heatmap_$type.pdf";
my $alldetected = "$outdir/alldetected_$type.tsv";
my $selected_log2fc = "$outdir/selected_log2fc_$type.tsv";
my $sessioninfo = "$outdir/sessionInfo.txt";
my $inputfile=$inputdir."/".$type.$input_file_suffix;
my $col=1;
$col=2 if ($type eq "isoforms");

my $heatmapR="";
if (lc($heatmap) eq "yes")
{
$heatmapR=qq/
  f1<-res[!is.na(res\$padj) & !is.na(res\$log2FoldChange), ]
  res_selected<-f1[(f1\$padj<$padj_cutoff & abs(f1\$log2FoldChange)>log2($foldChange_cutoff)),]
  if (dim(res_selected)>4 && length(conds)>2)
  {
  ld <- log(filtd[rownames(res_selected),]+0.1,base=2)

  cldt <- scale(t(ld), center=TRUE, scale=TRUE);
  cld <- t(cldt)

  dissimilarity <- 1 - cor(cld)
  distance <- as.dist(dissimilarity)
  pdf("$pdfname")
  plot(hclust(distance),  main="Dissimilarity = 1 - Correlation", xlab="")
  heatmap.2(cld, col=redgreen(75), scale="row",
          key=TRUE, symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(12,8),trace="none",srtCol=45)
  dev.off()
  }
/;
}

open(OUT, ">$output");
my $rscript=qq/
library("DESeq2")
library("ggplot2")
library("gplots")
analyseDE <-  function(data,cond, fitType, tablefile )
{
  tryCatch({
  data1<-data.frame(data)
  cols = c(1:dim(data1)[2]);
  data1[,cols] = apply(data1[,cols], 2, function(x) as.numeric(as.integer(x)))
  conds <- factor(cond)
  colData = as.data.frame((colnames(data1)));
  colData<-cbind(colData, conds)
  colnames(colData) = c("Cond","group");
  groups = factor(colData[,2]);
  sumd = apply(X=data1,MARGIN=1,FUN=sum);

  filtd = subset(data1, sumd > 10);

  dds = DESeqDataSetFromMatrix(countData=as.matrix(filtd), colData=colData, design = ~ group);

  dds <- DESeq(dds);

  res <- results(dds);
  
  write.csv(res, file=tablefile)    
  $heatmapR
  ult<-cbind(data[rownames(res_selected), ], res[rownames(res_selected), c("padj", "log2FoldChange")],  2 ^ (res[rownames(res_selected), "log2FoldChange"]) )
  colnames(ult)[dim(ult)[2]]<-"foldChange"
  write.table(ult, "$selected_log2fc", sep="\t", col.names=NA)

  all<-cbind(data[rownames(filtd), ], res[rownames(res), c("padj", "log2FoldChange")], 2 ^ (res[rownames(res), "log2FoldChange"]) )
  colnames(all)[dim(all)[2]]<-"foldChange"
  write.table(all, "$alldetected", sep="\t", col.names=NA)
  },
  finally={
     sessionInfo()
  })
}

file<-"$inputfile"
dat<- data.frame(read.table(file,sep="\\t", header=TRUE, row.names=$col, quote = "\\"", dec = "."), stringsAsFactors=TRUE);

data <- dat[, $cols]
cond <- factor( $conds )
analyseDE(data, cond, "$fitType","$table")
/;
print $rscript; 


print OUT $rscript; 
close(OUT);

my $com="$rscriptCMD $output > $sessioninfo 2>&1";
`$com`;
#die "Error 22: Cannot run Rscript:" if ($?);
my $verstring =`grep DESeq2_ $sessioninfo`;
$verstring =~/(DESeq[^\s]+)/;
my $deseq_ver=$1;
$deseq_ver = "DESeq2" if ($deseq_ver =~/^$/);
   
$com="sed -i 's/\"\"/name/' $selected_log2fc 2>/dev/null && sed -i 's/\"\"/name/' $alldetected 2>/dev/null &&";
$com.="sed -i 's/\"//g' $selected_log2fc 2>/dev/null && sed -i 's/\"//g' $alldetected 2>/dev/null && ";
$com.="echo \"$wkey\t$deseq_ver\tdeseq\t$deseqdir/alldetected_$type.tsv\" >> $puboutdir/reports.tsv && ";
$com.="echo \"$wkey\t$deseq_ver\tdeseq\t$deseqdir/selected_log2fc_$type.tsv\" >> $puboutdir/reports.tsv && ";
$com.="echo \"$wkey\t$deseq_ver\tdeseq\t$deseqdir/rscript_$type.R\" >> $puboutdir/reports.tsv ";
if (lc($heatmap) eq "yes")
{
  $com.="&& echo \"$wkey\t$deseq_ver\tdeseq\t$deseqdir/heatmap_$type.pdf\" >> $puboutdir/reports.tsv ";
}
`$com`;
#die "Error 21: Cannot run DESeq2 output files:" if ($?);
}

sub checkCols
{
 my $cols=$_[0];
 my @arr = split(/,/, $cols);
 print $cols."\n";
 foreach my $val (@arr)
 {
    if ($val=~/^[\d]+/){
      $cols=~s/$val/X$val/g;
    }
 }
  
 return $cols;
}

__END__


=head1 NAME

stepDESeq2.pl

=head1 SYNOPSIS  

stepDESeq2.pl -i input <fastq> 
            -o outdir <output directory> 

stepDESeq2.pl -help

stepDESeq2.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -i  input file <fastq format> 
    
=head2 -o outdir <output directory>

the output files will be "$outdir/" 


=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program runs DESeq2 in R

=head1 EXAMPLE


stepDESeq2.pl -c cols
            -o ~/out

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


