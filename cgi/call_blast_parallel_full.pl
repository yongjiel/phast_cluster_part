#!/usr/bin/perl -w
#usage : perl call_blast_parallel_full.pl <full_path_faa_file> <full_path_DB>

# call blast in pararllel
use File::Basename;

my $core_number = 56;
my $time1= time;
my $cgi_dir = "/home/prion/phage/cgi";
my $faa_file = $ARGV[0];
my $log_file = "$faa_file.log";
if (-e $log_file){
	system("rm -rf $log_file");
}
my $dir = dirname($faa_file);
my $faa_file_basename= basename($faa_file);

my $database = $ARGV[1];

# cut faa file into pieces.
my @arr = split (">", `cat $faa_file`);
my $piece=0;
for(my $i=0; $i <=$#arr; $i++){
  next if ($arr[$i] eq '' or $arr[$i]=~/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/s);
  my $sufix  = ($i % $core_number)+1;
  open(OUT, ">> $faa_file\_$sufix");
  print OUT '>'. $arr[$i];
  close OUT;
  if ($sufix >$piece){
    $piece=$sufix;
  }
}
if ($piece ==0){
  system("echo 'Nothing in .faa file. Program exit!' >>$log_file");
  exit(0);
}else{
	system("echo 'Split .faa file into $piece!' >> $log_file");
}

# call single_node_blast
chdir $dir;
system("echo change to $dir >>$log_file");
system("echo 'qsub -t 1-$piece  -q all.q -sync yes  $cgi_dir/full_single_blast.pl  $faa_file_basename $database' >> $log_file");
system("qsub -t 1-$piece  -q all.q -sync yes  $cgi_dir/full_single_blast.pl  $faa_file_basename $database")==0 or system("echo $! >> $log_file");; 

#combine all the files
my $file_str='';
for (my $i = 1; $i<=$piece; $i++){
	$file_str .= "$faa_file_basename\_$i\_out  ";
}
system("echo '$file_str' >>$log_file");
system("cat $file_str > $faa_file_basename.blast_out");
my $time2=time;
my $time_diff=$time2-$time1;
system("echo 'time_diff=$time_diff' >> $log_file");
system("echo 'program exit' >> $log_file");
exit;


		
