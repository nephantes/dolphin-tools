#!/usr/bin/env perl
#########################################################################################
#                                       rerun.pl
#########################################################################################
#        This program removes step-files to indicate from what step to rerun
#########################################################################################
# AUTHORS:
# Hennady Shulha, PhD 
# Alper Kucukural, PhD 
#########################################################################################
##################################### LIBRARIES AND PRAGMAS #############################
use Getopt::Long;
######################################## PARAMETER PARSING ##############################
GetOptions( 
	'outdir=s'            => \$outdir,
	'rerun=s'           => \$rerun,
	'servicename=s'       => \$servicename,
	'jobsubmit=s'    => \$jobsubmit,
) or die("Deleter step got unrecognized options.\n");
########################################## MAIN PROGRAM ##################################
if($rerun eq "0"){}
elsif($rerun eq "1"){`rm $outdir/tmp/track/Bam2BW*;rm $outdir/tmp/track/igvConverter*;rm $outdir/tmp/track/StepMapping*;rm $outdir/tmp/track/stepMerge*;rm $outdir/tmp/track/stepSplitFastQ*;rm $outdir/tmp/track/submitAgg*;rm $outdir/tmp/track/submitMacs*;`}
elsif($rerun eq "2"){`rm $outdir/tmp/track/Bam2BW*;rm $outdir/tmp/track/igvConverter*;rm $outdir/tmp/track/submitAgg*;rm $outdir/tmp/track/submitMacs*;`}
elsif($rerun eq "3"){`rm $outdir/tmp/track/Bam2BW*;rm $outdir/tmp/track/submitAgg*;rm $outdir/tmp/track/submitMacs*;`}
elsif($rerun eq "4"){`rm $outdir/tmp/track/submitAgg*;rm $outdir/tmp/track/submitMacs*;`}
elsif($rerun eq "5"){`rm $outdir/tmp/track/submitAgg*;`}
else
{
 print STDERR "Weired input in Rerun field\n";
 exit(1);
}

__END__

=head1 NAME

rerun.pl

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
