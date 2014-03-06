#!/usr/bin/env perl
#########################################################################################
#                                       mapping.pl
#########################################################################################
#                    This program does Bowtie1 mapping, conversion into BAM
#########################################################################################
# AUTHORS:
# Hennady Shulha, PhD 
# Alper Kucukural, PhD 
#########################################################################################
####################################### LIBRARIES AND PRAGMAS ###########################
use Getopt::Long;
####################################### PARAMETER PARSUING ##############################
GetOptions( 
	'outdir=s'            => \$outdir,
	'input=s'             => \$inputdir,
	'program=s'           => \$param,
	'bowtie=s'           => \$bowtie,
	'mapcommand=s'	=> \$mapper,
	'formatcommand=s'	=> \$formatter,
	'geneindex=s'           => \$geneindex,
	'servicename=s'       => \$servicename,
	'jobsubmit=s'    => \$jobsubmit,
) or die("Fastq_split step got unrecognized options.\n");
######################################### MAIN PROGRAM ##################################
if (! -e "$outdir")
{
 $res=`mkdir -p $outdir`;
 if($res != 0)
 {
  print STDERR "Failed to create output folder for Bowtie mapped files\n";
  exit(1);
 }
}
opendir $dir, "$inputdir" or die "Cannot open directory with files for mapping!";
@files = readdir $dir;
closedir $dir;
foreach $e(@files)
{
  next if ($e =~ m/^\./);
  $fastq="$inputdir/$e";
  if($bowtie eq "bowtie")
  {
   $com="$mapper $param $geneindex $fastq | awk \'\{if(\\\$4 != 0) print(\\\$0)\}\' > $outdir/$e.sam\; $formatter view -b -S $outdir/$e.sam > $outdir/$e.bam\; $formatter sort $outdir/$e.bam $outdir/$e.sorted";
   $job=$jobsubmit." -n ".$servicename."_".$e." -c \"$com\"";
   $res=`$job`;
   print $job."\n";
   if($res != 0)
   {
    print STDERR "Failed to submit mapping for $e\n";
    exit(1);
   }   
  }
}

__END__

=head1 NAME

mapping.pl

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
