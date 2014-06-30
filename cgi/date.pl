#!/usr/bin/perl -w
#$ -S /usr/bin/perl
#$ -cwd
#$ -r yes
#$ -j no
#$ -pe smp 1

my $date = `date`;
my $id = $ENV{SGE_TASK_ID}; # Batch-scheduler assigns this.
system("echo $date > date\_$id.log");
exit;

