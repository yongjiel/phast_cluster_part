#!/usr/bin/perl -w
# this program uses one.q and calls sing_blast.pl in qsub
# call blast in pararllel
use File::Basename;

my $core_number = 56;
my $cgi_dir = "/home/prion/phage/cgi";
my $faa_file = $ARGV[0];
my $log_file = "$faa_file.log";
if (-e $log_file){
	system("rm -rf $log_file");
}
my $dir = dirname($faa_file);
my $faa_file_basename= basename($faa_file);

my $database = $ARGV[1];

my $last_node = 10; # max cluster node number to use

# load optional arguments used for testing
if (defined($ARGV[2])) {
	$last_node = $ARGV[2];
}
if (defined($ARGV[3])) {
	$core_number = $ARGV[3];
}

# cut faa file into pieces.
my @arr = split (">", `cat $faa_file`);
@arr = grep ($_ !~ /^[\n\s]*$/, @arr);
my $piece=0;
for(my $i=0; $i <=$#arr; $i++){
	my $sufix  = ($i % $core_number)+1;
	open(OUT, ">> $faa_file\_$sufix");
	print OUT '>'. $arr[$i];
	close OUT;
	if ($sufix >$piece){
		$piece=$sufix;
	}
}
if ($piece == 0){
	system("echo 'Nothing in .faa file. Program exit!' >>$log_file");
	exit(0);
}else{
	system("echo 'Split .faa file into $piece!' >> $log_file");
}

# Generate string specifying which nodes to run on.
my $node_set = '';
for (my $i = 1; $i <= $last_node; $i++) { # restrict to particular nodes
	$node_set .= 'botha-w' . $i;
	if ($i != $last_node) {
		$node_set .= '|';
	}
}

# call single_node_blast
chdir $dir;
system("echo change to $dir >>$log_file");
system("echo 'qsub -t 1-$piece  -pe smp 1 -q one.q -l h=\"$node_set\" -sync yes  $cgi_dir/single_blast_old_PHAST_for_testing.pl  $faa_file_basename $database' >> $log_file");
system("qsub -t 1-$piece  -pe smp 1 -q one.q -l h=\"$node_set\" -sync yes  $cgi_dir/single_blast_old_PHAST_for_testing.pl  $faa_file_basename $database")==0 or system("echo $! >> $log_file");; 

#combine all the files
my $file_str='';
for (my $i = 1; $i<=$piece; $i++){
	$file_str .= "$faa_file_basename\_$i\_out  ";
}
system("echo '$file_str' >>$log_file");
system("cat $file_str > $faa_file_basename\_blast_out");
system("echo 'program exit' >> $log_file");
exit;


		
