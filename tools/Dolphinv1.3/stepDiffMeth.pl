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
 my $stepsize         = "";
 my $bedfile          = "";
 my $topN             = "";
 my $mincoverage      = "";
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
    'outdir=s'       => \$outdir,
    'name=s'         => \$name,
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
my $inputdir   = "$outdir/methylKit";
$outdir   = "$outdir/diffmeth_".$name;
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

runMethylKit($inputdir, $bedfile, $samplenames, $conds, $outdir, $puboutdir, $wkey);

`cp -R $outdir $puboutdir/.`;

sub runMethylKit
{
my ($inputdir, $bedfile, $samplenames, $conds, $outdir, $puboutdir, $wkey)=@_;
my $output = "$outdir/rscript_$name.R";
my $sessioninfo = "$outdir/sessionInfo.txt";

open(OUT, ">$output");
my $rscript = qq/
library("methylKit")
  inputdir<-$inputdir; samplenames<-$samplenames; conds<-$conds; outdir<-"$outdir";
  conds<-conds-1
  bedfile<-"$bedfile"
  
  load(file=paste(inputdir,"\/calcdata.rda"))
  
  tiles_comp=reorganize(meth_tiles,sample.ids=samplenames,
                        treatment=conds )
  
  myDiff<-calculateDiffMeth(tiles_comp,slim=TRUE,weigthed.mean=TRUE,num.cores=4)
  
  myDiff_25p.hyper<-get.methylDiff(myDiff,difference=1,qvalue=0.01,type="hyper")
  myDiff_25p.hypo<-get.methylDiff(myDiff,difference=1,qvalue=0.01,type="hypo")
  
  data<-getData(meth_tiles)
  data_cpgcov<-getData(meth_tiles_cpgcov)
  
  difftiles.hyper<-getData(myDiff_25p.hyper)
  difftiles.hypo<-getData(myDiff_25p.hypo)
  difftiles<-getData(myDiff)
  
  write.table(difftiles.hyper, paste(outdir,"\/difftiles.hyper.tsv",sep=""))
  write.table(difftiles.hypo, paste(outdir,"\/difftiles.hypo.tsv",sep=""))
  write.table(difftiles, paste(outdir,"\/difftiles.tsv",sep=""))
  
  gene.obj=read.transcript.features(bedfile)
  ann<-annotate.WithGenicParts(myDiff,gene.obj)
  
  pdf(paste(outdir, "\/geneannot.pdf", sep=""))
  plotTargetAnnotation(ann,precedence=TRUE)
  dev.off()
  
  promoters=regionCounts(meth_tiles,gene.obj\$promoters)
  write.table(promoters, paste(outdir,"\/promoters.tsv",sep=""))
/;
print $rscript; 


print OUT $rscript; 
close(OUT);

my $com="$rscriptCMD $output > $sessioninfo 2>&1";

my $job=$jobsubmit." -n ".$servicename."_".$name." -c \"$com\"";
print $job."\n";   
`$job`;
die "Error 25: Cannot run the job:".$job if ($?);
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


