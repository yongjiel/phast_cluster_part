#!/usr/bin/perl -w
use  Cwd;

my $exec_dir = "/home/prion/phage/cgi";
my $tmp_dir = "/home/prion/phage/DB/temp_bac";

my $start_time = time;
chdir $tmp_dir;
print  `date`;
system("rm -rf *") if (-e "list");

system("wget 'ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.faa.tar.gz'  2>&1 |cat >/dev/null");
if (-e "all.faa.tar.gz"){
	print "all.faa.tar.gz is downloaded!\n";
}else{
	print "No all.faa.tar.gz is downloaded!\n";
	exit;
}
system("tar -xvzf all.faa.tar.gz 2>&1");
system("ls >list");

open(IN, "list") or die "Cannot open list";
my $output="bacteria_all.db";
unlink ($output);
open(OUT, ">$output") or die "Cannot write $output";
my $f =0;
my $last ='';
my $curr_dir =getcwd;
my $count=0;
while(<IN>) {
	chomp ($_);
	if ($_ eq 'list' or  $_ eq 'make_bac_database.pl' or $_ eq 'all.faa.tar.gz' or  $_ =~/$output/  ){
		next;
	}
	
	chdir $_;
	my $v_name = $_;
	
	my $file_list =`ls`;
	my @arr=split("\n", $file_list);
	foreach my $file (@arr){ 
		chomp($file); 
		$count++;
		#print "     $v_name/$file\n";
		my $data = `cat $file`;
		print OUT $data;
		$f=1;
	}
	chdir $curr_dir;

}
close IN;
if ($f==1){
	print "$output rewritten\n";
}
close OUT;
print "In all.faa.tar.gz, we found $count phage organisms!\n";


system("grep '>' $output > $output\_list");
if (-s "$output"){
	#print "Copy Bacteria_all* to upper level\n";
	#system("cp Bacteria_all* ../.");
	# propogate to prion\@botha1.cs.ualberta.ca
	#print "Propogate Bacteria_all* to prion\@botha1.cs.ualberta.ca\n";
	#system("$exec_dir/propogate_database.sh  /home/phast/public_html/phage_finder/DB/temp_bac/$output");
	#cleanup
	#system("rm -rf *") if (-s "list");
}else{
	print "$output  has no content. Please check!\n";
}
my $end_time = time;
my $time = $end_time - $start_time;
my $min = int($time /60);
my $sec = $time % 60;
print "running time = $min min $sec sec\n";
print "Program exit!\n\n\n";
exit;


