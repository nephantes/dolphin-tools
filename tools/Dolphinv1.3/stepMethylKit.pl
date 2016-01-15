#!/usr/bin/env perl

#########################################################################################
#                                       stepMethylKit.pl
#########################################################################################
# 
#  This program  run methylKit in R.
#
#
#########################################################################################
# AUTHORS:
#
# Alper Kucukural, PhD 
# Jan 14, 2016
#########################################################################################


############## LIBRARIES AND PRAGMAS ################

 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 

#################### VARIABLES ######################
 my $samplenames      = "";
 my $conds            = "";
 my $gbuild           = "";
 my $outdir           = "";
 my $strand           = "";
 my $tilesize         = "";
 my $bedfile          = "";
 my $rscriptCMD       = "";
 my $pubdir           = "";
 my $wkey             = "";
 my $name             = "";
 my $jobsubmit        = "";
 my $servicename      = "";
 my $help             = "";
 my $print_version    = "";
 my $version          = "1.0.0";
################### PARAMETER PARSING ####################

my $cmd=$0." ".join(" ",@ARGV); ####command line copy

GetOptions( 
    'samplenames=s'  => \$samplenames,
    'conds=s'        => \$conds,
    'gbuild=s'       => \$gbuild,
    'outdir=s'       => \$outdir,
    'strand=s'       => \$strand,
    'name=s'         => \$name,
    'tilesize=s'     => \$tilesize,
    'bedfile=s'      => \$bedfile,
    'rscriptCMD=s'   => \$rscriptCMD,
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

pod2usage( {'-verbose' => 0, '-exitval' => 1,} ) if ( ($samplenames eq "") or ($conds eq "") or ($outdir eq "") );	

################### MAIN PROGRAM ####################
#    maps the reads to the ribosome and put the files under $outdir/after_ribosome directory

my $inputdir = "$outdir/mcall";
my $input_file_suffix = ".methylkit.txt";


$outdir   = "$outdir/meth_".$name;
`mkdir -p $outdir`;
$samplenames=~s/[\s\t]+//g;
$conds=~s/[\s\t]+//g;
$samplenames=~s/,/\",\"/g;
$samplenames="c(\"$samplenames\")";
$conds=~s/,/\",\"/g;
$conds=~s/Cond//g;
$conds=~s/\"//g;
$conds="c($conds)";

if (lc($strand) !~/^no/) {
    $strand = "T";
}
else{
    $strand ="F";
}


my $puboutdir   = "$pubdir/$wkey";
`mkdir -p $puboutdir`;

runMethylKit($inputdir, $input_file_suffix, $samplenames, $conds, $gbuild, $outdir, $strand, $tilesize, $puboutdir, $wkey);

`cp -R $outdir $puboutdir/.`;

sub runMethylKit
{
my ($inputdir, $input_file_suffix, $samplenames, $conds, $gbuild, $outdir, $strand, $tilesize, $puboutdir, $wkey)=@_;
my $output = "$outdir/rscript_$name.R";
my $sessioninfo = "$outdir/sessionInfo.txt";
print "inputdir<-\"$inputdir\";input_file_suffix<-\"$input_file_suffix\"; samplenames<-$samplenames; conds<-$conds; gbuild<-\"$gbuild\"; outdir<-\"$outdir\"; strand<-$strand; tilesize<-$tilesize;";

open(OUT, ">$output");
my $rscript = qq/
push <- function(l, ...) c(l, list(...))
outerJoin <- function(data1, data2,data3, fields)
{
  d1 <- merge(data1, data2, by=fields, all=TRUE)
  d2 <- merge(data3, d1,  by=fields, all=TRUE)
  d2[,fields]
}
runMethylSeq <- function(inputdir, input_file_suffix, samplenames, conds, gbuild, outdir, strand, tilesize)
{
  statspdf<-paste(outdir,"\/stats.pdf", sep="")
  analysispdf<-paste(outdir,"\/analysis.pdf", sep="")
  file.list<-list()

  for (i in seq(samplenames))
  {
    file.list<-push(file.list, paste(inputdir, "\/", samplenames[i], input_file_suffix, sep=""))
  }
  snames<-list()
  for (i in seq(samplenames))
  {
    snames<-push(snames, samplenames[i])
  }
  
  myobj=read( file.list,
              sample.id=snames,assembly=gbuild,treatment=conds)
  
  myobj.cpgcov<-myobj
  for (i in seq(samplenames))
  {
    myobj.cpgcov[[i]]\$coverage<-1
  }
  
  pdf(statspdf) 
  for (i in seq(samplenames))
  {
    getMethylationStats(myobj[[i]],plot=T,both.strands=F)
    getCoverageStats(myobj[[i]],plot=T,both.strands=F)
  }
  
  dev.off()
  meth=unite(myobj)
  
  tiles<-tileMethylCounts(myobj,win.size=tilesize,step.size=tilesize)
  tiles_cpgcov<-tileMethylCounts(myobj.cpgcov,win.size=tilesize,step.size=tilesize)
  
  meth_tiles<-unite(tiles)
  meth_tiles_cpgcov<-unite(tiles_cpgcov)
  
  pdf(analysispdf) 
  getCorrelation(meth_tiles,plot=T)
  clusterSamples(meth_tiles, dist="correlation", method="ward", plot=T)
  PCASamples(meth_tiles, screeplot=T)
  PCASamples(meth_tiles,screeplot=FALSE, adj.lim=c(0.0004,0.1),
             scale=TRUE,center=TRUE,comp=c(1,2),transpose=TRUE,sd.filter=TRUE,
             sd.threshold=0.5,filterByQuantile=TRUE,obj.return=FALSE)
  dev.off()
  
  tiles_comp=reorganize(meth_tiles,sample.ids=samplenames,
                        treatment=conds )
  
  myDiff<-calculateDiffMeth(tiles_comp,slim=TRUE,weighted.mean=TRUE,num.cores=4)
  
  myDiff_25p.hyper<-get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
  myDiff_25p.hypo<-get.methylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
  
  data<-getData(meth_tiles)
  data_cpgcov<-getData(meth_tiles_cpgcov)
  
  difftiles.hyper<-getData(myDiff_HF_LP_25p.hyper)
  difftiles.hypo<-getData(myDiff_HF_LP_25p.hypo)
  difftiles<-getData(myDiff_HF_LP_25p)
  
  write.table(difftiles.hyper, paste(outdir,"\/difftiles.hyper.tsv",sep=""))
  write.table(difftiles.hypo, paste(outdir,"\/difftiles.hypo.tsv",sep=""))
  write.table(difftiles, paste(outdir,"\/difftiles.tsv",sep=""))
  
  rownames(data)<-paste(data\$chr,data\$start, data\$end,sep="_")
  
  cols<-c()
  for (i in seq(samplenames))
  {
    cols<-c(cols, paste("coverage", i, sep="" ) )
  }
  
  norm_data<-cbind(rowSums(data[, cols]),
                   rowSums(data_cpgcov[, cols]))
  for (i in seq(samplenames))
  {
    norm_data<-cbind(norm_data, data[,paste("numCs", i, sep="" )]\/data[, paste("coverage", i, sep="" ) ] )
  }
  
  rownames(norm_data) <- rownames(data)
  colnames(norm_data) <- c("Cov", "Cpg", samplenames)
  snames<-samplenames
  
  filtmaxcoverage<-cbind(apply(norm_data[norm_data[,"Cov"]\/norm_data[,"Cpg"]>5,3:dim(norm_data)[2]], 1, function(x) max(x)),1)
  
  plot(density( filtmaxcoverage[,1]))
  lowelim<-norm_data[filtmaxcoverage[,1]>0.6, ]
  
  write.table(lowelim, paste(outdir,"\/after_elimination.tsv",sep=""))
  
  cv<-cbind(apply(lowelim, 1, function(x) (sd(x,na.rm=TRUE)\/mean(x,na.rm=TRUE))), 1)
  
  #cv<-cbind(apply(norm_data[norm_data[,"Cov"]\/norm_data[,"Cpg"]>5 ,snames], 1, function(x) (sd(x,na.rm=TRUE))), 1)
  colnames(cv)<-c("coeff", "a")
  
  withcvlowelim<-cbind(cv[,1], lowelim)
  colnames(withcvlowelim)[1]<-"Coeff"
  write.table(withcvlowelim, paste(outdir,"\/after_elimination_with_coeff.tsv",sep=""))
  
  cvsort<-cv[order(cv[,1],decreasing=TRUE),]
  
  cvsort_top2000 <- cvsort[1:2000,]
  
  selected<-data.frame(norm_data[rownames(cvsort_top2000),snames])
  #write.table(combData[rownames(cvsort_top2000),], paste(outdir,"\/combData_cov_top2000.tsv",sep=""))
  colnames(selected) <- snames
  write.table(selected, paste(outdir,"\/most_cv_top2000.tsv",sep=""))
  gene.obj=read.transcript.features(bedfile)
  ann<-annotate.WithGenicParts(myDiff,gene.obj)
  
  pdf(paste(outdir, "\/geneannot.pdf", sep=""))
  plotTargetAnnotation(ann,precedence=TRUE)
  dev.off()
  
  promoters=regionCounts(meth_tiles,gene.obj\$promoters)
  write.table(promoters, paste(outdir,"\/promoters.tsv",sep=""))
}

runMethylSeq <- function("$inputdir","$input_file_suffix", $samplenames, $conds, "$gbuild", "$outdir", $strand, $tilesize)
/;
print $rscript; 


print OUT $rscript; 
close(OUT);

my $com="$rscriptCMD $output > $sessioninfo 2>&1";
`$com`;
die "Error 22: Cannot run Rscript:" if ($?);
my $verstring =`grep methylKit_ $sessioninfo`;
$verstring =~/(methylKit[^\s]+)/;
my $methtylit_ver=$1;
#$com.="echo \"$wkey\t$methtylit_ver\tdeseq\t$deseqdir/alldetected_$type.tsv\" >> $puboutdir/reports.tsv && ";
#$com.="echo \"$wkey\t$methtylit_ver\tdeseq\t$deseqdir/selected_log2fc_$type.tsv\" >> $puboutdir/reports.tsv && ";
#$com.="echo \"$wkey\t$methtylit_ver\tdeseq\t$deseqdir/rscript_$type.R\" >> $puboutdir/reports.tsv ";
#`$com`;
#die "Error 21: Cannot run DESeq2 output files:" if ($?);
}

__END__


=head1 NAME

stepMethylKit.pl

=head1 SYNOPSIS  

stepMethylKit.pl -sa samplenames <comma separated> 
            -o outdir <output directory> 

stepMethylKit.pl -help

stepMethylKit.pl -version

For help, run this script with -help option.

=head1 OPTIONS

=head2 -sa  samplnames 

{samples}.mehtylkit.txt will be input files.
    
=head2 -o outdir <output directory>

the output files will be "$outdir" 


=head2 -help

Display this documentation.

=head2 -version

Display the version

=head1 DESCRIPTION

 This program runs MethylKit in R

=head1 EXAMPLE


stepMethylKit.pl -c conds
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


