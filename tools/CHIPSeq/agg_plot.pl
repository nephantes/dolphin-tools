#!/usr/bin/env perl
#########################################################################################
#                                       agg_plot.pl
#########################################################################################
#                    This program perform ACT aggregation for the mapped tags
#########################################################################################
# AUTHORS:
# Hennady Shulha, PhD 
# Alper Kucukural, PhD 
#########################################################################################
####################################### LIBRARIES AND PRAGMAS ###########################
use Getopt::Long;
use File::Basename;
####################################### PARAMETER PARSUING ##############################
GetOptions(
	'program=s'         => \$program,
	'genome=s'       => \$genome, #(hg19.chromInfo)
	'reference=s'    => \$reference, #(refseq4col)
 	'outdir=s'       => \$indir,
 	'bedtoolsgencov=s'       => \$bedtoolsgencov,
 	'creation=s'       => \$creationpdf,
	'input=s'        => \$input,
	'jobsubmit=s'    => \$jobsubmit,
	'servicename=s'  => \$servicename,
) or die("Unrecognized optioins for ACT aggregator.\n");
################### MAIN PROGRAM ####################
$outdir="$indir/aggregationout";
if (! -e "$outdir")
{
 `mkdir -p $outdir`;
 if($? != 0)
 {
  print STDERR "Failed to create output folder for aggregation plot output\n";
  exit(1);
 }
}
@prefiles=split(/:/,$input);
for($i=0;$i<@prefiles;$i++) 
{
 @file=split(/,/,$prefiles[$i]);
 if (@file eq 2)
 {
  ($filename)  = fileparse($file[0]);
  $com="$bedtoolsgencov -bga -ibam $indir/$filename.sorted.bam -g $genome > $outdir/$filename.bed;\n";
  $com.="awk '{print \\\$1\\\"\\\\t\\\"\\\$2\\\"\\\\t\\\"\\\$4}' $outdir/$filename.bed > $outdir/$filename.sig;\n";
  $com.="python $program --nbins=50 --mbins=50 --radius=5000 --region $reference $outdir/$filename.sig > $outdir/$filename.aggregation_plot.out;\n";
  $com.="$creationpdf --args $outdir/$filename.aggregation_plot.out;\n";
  $a="$filename.agg";
  $job=$jobsubmit." -n ".$servicename."_".$a." -c \"$com\"";
  `$job`;
  print $job."\n";
  if($? != 0)
  {
   print STDERR "Failed to submit aggregation plot program\n";
   exit(1);
  } 
 }
}

__END__

=head1 NAME

agg_plot.pl

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