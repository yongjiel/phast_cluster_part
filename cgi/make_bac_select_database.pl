#!/usr/bin/perl -w
use  Cwd;
use File::Copy;
use File::Path;

my $exec_dir = "/home/prion/phage/cgi";
my $tmp_dir = "/home/prion/phage/DB/temp_bac_select";
my $tmp_genome_dir='/home/prion/phage/DB/temp_genome';
if (-d $tmp_genome_dir){
	chdir $tmp_genome_dir;
	system("rm -rf *");
	system("wget 'ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.gbk.tar.gz' 2>&1 |cat >/dev/null");
	system("tar -xvzf all.gbk.tar.gz");
}
if (!-d $tmp_dir){
	print "There is no $tmp_dir, program exit!\n";
	exit(-1);
}
chdir $tmp_dir;

my $start_time = time;
print  `date`;
system("rm -rf *") if (-e "list");

system("wget 'ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.faa.tar.gz' 2>&1 |cat >/dev/null");
if (-e "all.faa.tar.gz"){
	print "all.faa.tar.gz is downloaded!\n";
}else{
	print "No all.faa.tar.gz is downloaded!\n";
	exit(-1);
}

system("tar -xvzf all.faa.tar.gz > zz_tar_list_file");
system("perl $exec_dir/get_genome_NC_protein_gis.pl  $tmp_dir"); #create z_genome_NC_protein_gis file;
copy  "z_genome_NC_protein_gis", "../z_genome_NC_protein_gis";
system("perl $exec_dir/get_organism_NC_match_NC_only.pl"); # create z_current_organism_NC_list and z_current_NC_list file in parent directory.
system("ls >list");

open(IN, "list") or die "Cannot open list";
my $output="bacteria_all_select.db";
unlink ($output);
open(OUT, ">$output") or die "Cannot write $output";
my $last ='';
my $curr_dir =getcwd;
my $count=0;
while(<IN>) {
	chomp ($_);
	if ($_ eq 'list' or  $_ eq 'make_bac_select.pl' or $_ eq 'all.faa.tar.gz' or  $_ =~/bacteria_all_select\.db/  ){
		next;
	}
	
	chdir $_;
	my $v_name = $_;
	
	my $file_list =`ls`;
	my @arr=split("\n", $file_list);
	my $file=$arr[0];
		chomp($file); 
		$count++;
		#print "     $v_name/$file\n";
		my $data = `cat $file`;
		print OUT $data;
	
	chdir $curr_dir;

}
close IN;
close OUT;
print "In all.faa.tar.gz, we found $count bacteria organisms!\n";

open(IN, "bacteria_all_select.db") or die "Cannot open bacteria_all_select.db";
open(OUT, ">bacteria_all_select.db.tmp") or die "Cannot write bacteria_all_select.db.tmp"; 
my %hash =();
my $flag = 0;
while(<IN>) {
	if ($_=~/>gi\|(\d+)/){
		#>gi|158333234|ref|YP_001514406.1| NUDIX hydrolase [Aca
		$hash{$1} += 1;
		if ($hash{$1}>1){
			$flag = 1;
		}else{
			$flag =0;
			print OUT $_;
		}
	}else{
		if ($flag==0){
			print OUT $_;
		}
	}
}
close IN;
close OUT;
system("mv -f bacteria_all_select.db.tmp  bacteria_all_select.db");
my $error_flag = 0;
if (-s "bacteria_all_select.db"){
	system("echo 'make bacteria_all_select_header_lines.db'");
	system("grep '>' bacteria_all_select.db > bacteria_all_select_header_lines.db");
	system("rm -rf bacteria_all_select.db.*");
	system(" /home/prion/blast/bin/formatdb  -i  bacteria_all_select.db -o T ");
	print "Copy bacteria_all_select* to upper level\n";
	system("cp bacteria_all_select* ../.");
	#clean up
	system("rm -rf * ") if (-s 'list');
}else{
	print "bacteria_all_select.db  has no content. Please check!\n";
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


