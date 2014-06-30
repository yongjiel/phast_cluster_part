#/usr/bin/perl 
#
#$ -S /usr/bin/perl
#$ -cwd
#$ -r yes
#$ -j no
#$ -pe smp 1     
#
##===============================================================================
##
##         FILE:  date.pl
##
##        USAGE:  ./date.pl  
##
##  DESCRIPTION:  
##
##      OPTIONS:  ---
## REQUIREMENTS:  ---
##         BUGS:  ---
##        NOTES:  ---
##       AUTHOR:  YOUR NAME (), 
##      COMPANY:  
##      VERSION:  1.0
##      CREATED:  06/14/2011 11:01:50 AM 
##     REVISION:  ---
##===============================================================================
#
use strict; 
use warnings;

my $fna_file_basename = $ARGV[0];

my $id = $ENV{SGE_TASK_ID}; # Batch-scheduler assigns this.

my $glimmer_exec = "/home/prion/phage/glimmer3.02/scripts/g3-iterated.csh";
system("csh $glimmer_exec  $fna_file_basename\_$id   $fna_file_basename\_$id");
exit;
