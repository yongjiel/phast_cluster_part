#!/usr/bin/perl -w

# update database automatically
my $exec_dir ="/home/prion/phage/cgi";
my $DB_dir ="/home/prion/phage/DB";
my $node_num = 11;
my $Dest_DB_dir = "/usr/scratch/phage/DB";
my $ok1 = system("perl $exec_dir/make_virus_database.pl 2>&1|cat >> $DB_dir/make_virus_database.log");
my $ok2 = system("perl $exec_dir/make_prophage_virus_db.pl  2>&1|cat >> $DB_dir/make_prophage_virus_db.log");
my $ok3 = system("perl $exec_dir/make_bac_select_database.pl 2>&1|cat  >> $DB_dir/make_bac_select_database.log");
if ($ok1==0 && $ok2==0 && $ok3==0){
	my $date = `date`;
	my $update_logfile = "$DB_dir/update.log" ;
	system("echo '$date' > $update_logfile");
}
`perl $exec_dir/index_db.pl $DB_dir/bacteria_all_select_header_lines.db`;
# copy files to destination dir in each node.
for (my $i = 1; $i<= $node_num; $i++){
	`scp -i ~/.ssh/scp-key $DB_dir/bacteria_all_select.db* $DB_dir/bacteria_all_select_header_lines.db $DB_dir/prophage_virus.db* $DB_dir/prophage_virus_header_lines.db prion\@botha-w$i\:$Dest_DB_dir/.`;
}
exit;

