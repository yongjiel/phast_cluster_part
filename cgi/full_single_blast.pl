#!/usr/bin/perl -w
#$ -S /usr/bin/perl
#$ -cwd
#$ -r yes
#$ -j no
#$ -pe smp 1 


# run blast on a single node of our cluster.

my $faa_file_basename = $ARGV[0];
my $database = $ARGV[1];

my $id = $ENV{SGE_TASK_ID}; # Batch-scheduler assigns this.

my $blast_exec = '/home/prion/blast/bin/blastall';
my $command = "$blast_exec  -p blastp -d $database -m 8 -e 0.0001 -i $faa_file_basename\_$id -v 10 -b 10 -o $faa_file_basename\_$id\_out  -a 1 -F F ";
system("echo '$command'");
system($command) ;
exit;
