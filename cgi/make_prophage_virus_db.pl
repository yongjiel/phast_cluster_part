#!/usr/bin/perl -w

# make a db combining prophage and virus
my $exec_dir = "/home/prion/phage/cgi";
my $db_dir = "/home/prion/phage/DB";

print `date`;
my $start_time = time;
if (!(-e "$db_dir/virus.db")){
	print "$db_dir/virus.db not exist!\n";
	exit(-1);
}
system("cat $db_dir/temp_pro/prophage.db  $db_dir/virus.db > $db_dir/prophage_virus.db");

chdir $db_dir;
my $error_flag = 0;
if (-s "prophage_virus.db"){
	print "Make header lines file\n";
	system("grep '>' prophage_virus.db > prophage_virus_header_lines.db");
	print "formatdbing prophage_virus.db\n";
	system(" /home/prion/blast/bin/formatdb  -i    prophage_virus.db -s T");
}else{
	print "prophage_virus.db  has no content. Please check!\n";
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
