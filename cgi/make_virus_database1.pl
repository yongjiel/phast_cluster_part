#!/usr/bin/perl -w
use  Cwd;

my $exec_dir = "/home/prion/phage/cgi";
my $tmp_dir = "/home/prion/phage/DB/temp_vir";
if (!-d $tmp_dir){
	print "There is no $tmp_dir, program exit!\n";
	exit(-1);
}
my $start_time = time;
chdir $tmp_dir;
print  `date`;
system("rm -rf *") if (-e "list");
system("wget -O phage.html 'http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&opt=Virus&sort=genome' 2>&1 |cat >/dev/null");
if (-e "phage.html"){
	print "phage.html is generated!\n";
}else{
	print "No phage.html is generated!\n";
}
system("wget 'ftp://ftp.ncbi.nih.gov/genomes/Viruses/all.faa.tar.gz' 2>&1 |cat >/dev/null");
if (-e "all.faa.tar.gz"){
	print "all.faa.tar.gz is generated!\n";
}else{
	print "No all.faa.tar.gz is generated!\n";
	exit(-1);
}
system("tar -xvzf all.faa.tar.gz > /dev/null");
system("ls >list");
if (!(-s "phage.html")){
	print "Something wrong when wget -O phage.html 'http://www.ncbi.nlm.nih.gov/genomes/genlist.cgi?taxid=10239&type=6&name=Phages'\n";
}
my $need_phage = `cat phage.html`;
my @need_arr = $need_phage =~/>(NC_\d+)</gs;
print "Needed phage organism number from Phage.html =".scalar(@need_arr)."\n";
my @curr_arr=();
open(IN, "list") or die "Cannot open list";
my $output="virus.db";
unlink ($output);
open(OUT, ">$output") or die "Cannot write $output";
my $last ='';
my $curr_dir =getcwd;
my $count=0;
while(<IN>) {
	chomp ($_);
	if ($_ eq 'list' or $_ eq 'phage.html' or  $_ eq 'make_virus_database.pl' or $_ eq 'all.faa.tar.gz' 
	    or  $_ =~/virus\.db/  or  $_ =~/\.gbk/  or  $_ =~/\.ptt/){
		next;
	}
	
	chdir $_;
	my $v_name = $_;
	$v_name =~s/_uid.*//;
	my $prefix = get_prefix($v_name);
	my $file_list =`ls`;
	my @arr=split("\n", $file_list);
	foreach my $file (@arr){ 
		chomp($file);
		my $found = 0;
		if(scalar(@need_arr) !=0){
			foreach my $l (@need_arr){
				if ($file =~/$l/){
					$found = 1;
					last;
				}
			}
			if ($found ==0){
				next;
			} 
		}
		$count++;
		print "     $v_name/$file\n";
		my $data = `cat $file`;
		$data =~s/\[.*?\]/\[$v_name\]/gs;
		$data =~s/>gi/>$prefix-gi/gs;
		print OUT $data;
		$file =~s/\.faa//;
		push @curr_arr, $file;
	}
	chdir $curr_dir;

}
close IN;
if(scalar(@need_arr) !=0){
	print "For phage.html, in all.faa.tar.gz, we found $count phage organisms!\n";
}else{
	print "In all.faa.tar.gz, we found $count phage organisms!\n";
}
my $cou1=0;
for(my $i =0; $i <=$#need_arr; $i++){
	my $found =0;
	for(my $j = 0; $j <=$#curr_arr; $j++){
		if ($need_arr[$i] eq $curr_arr[$j]){
			$found = 1;
			last;
		}
	}
	if ($found ==0){
#		print "   $need_arr[$i] from need_arr not found in curr_arr, call genbank to get the gbk file!\n";
		system("perl $exec_dir/get_gbk.pl  $need_arr[$i]  $tmp_dir");# create .gbk fil
		if (-s "$tmp_dir/$need_arr[$i].gbk"){
			print "$need_arr[$i] back from genbank\n";
		}
		system("perl $exec_dir/gbk2faa.pl $need_arr[$i].gbk"); #create .faa file
		my $data = `cat $need_arr[$i].faa`;
		$data =~s/,\s*complement\(\d+\.\.\d+\)//gs;
		$data =~s/,\s*\d+\.\.\d+//gs;
		
		my @lines =split("\n", $data);
		foreach my $l (@lines){
			if ($l =~/>gi.*\[(\S+).*\s+(\S+)\]/ or $l =~/>gi.*\[(\S+)\]/){
				my $substr='';
				if (defined $2 ){
					$substr= "PHAGE_".substr($1,0,6)."_$2";
				}else{
					$substr= "PHAGE_".substr($1,0,6);
				}
				$substr=~s/[\s-]/_/g;
				$l =~s/>gi/>$substr-gi/;
			}
			print OUT $l."\n";
		}

	}
}
exit;
my $error_flag = 0;
if (-s "virus.db"){
	print "Copy virus.db* to upper level\n";
	system("cp virus.db* ../.");
	#cleanup 
	system("rm -rf *") if (-s "list");
}else{
	print "virus.db  has no content. Please check!\n";
	$error_flag = 1;
}
my $end_time = time;
my $time = $end_time - $start_time;
my $min = int($time /60);
my $sec = $time % 60;
print "running time = $min min $sec sec\n";
print "Program exit!\n\n\n";
if ($error_flag ==1){
	exit(-1);
}
exit;

sub get_prefix{
	my $v_name= shift;
	$v_name =~s/_phage//;
	$v_name =~s/Phage_//;	
	if ($v_name =~/^(\w+)_/){
		my $str = $1;
		$str_sub = substr($str, 0, 6);
		$v_name =~s/^(\w+)_/$str_sub\_/;
	}else{
		#print "Wierd, $v_name\n";
	}
	$v_name = "PHAGE_$v_name";		
	
	return $v_name;

}
