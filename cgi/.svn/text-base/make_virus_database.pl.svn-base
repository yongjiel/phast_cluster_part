#!/usr/bin/perl -w
use  Cwd;
use  File::Copy;

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
my %hash_old_virus = get_hash_old_virus_db();
print "Get ".scalar(keys %hash_old_virus)." NC from ../virus.db\n";
open(IN, "list") or die "Cannot open list";
my $output="virus.db";
unlink ($output);
open(OUT, ">$output") or die "Cannot write $output";
my $last ='';
my $curr_dir =getcwd;
my $count=0;
while(<IN>) {
	chomp ($_);
	next if (! -d $_);	
	chdir $_;
	my $v_name = $_;
	$v_name =~s/_uid.*//;
	my $prefix = get_prefix($v_name);
	my $file_list =`ls`;
	my @arr=split("\n", $file_list);
	foreach my $file (@arr){ 
		chomp($file);
		my $basename= $file;
		$basename=~s/\.faa//;
		next if (defined $hash_old_virus{$basename}); #means it is in old virus db
		for (my $i = 1; $i<=3; $i++){
			system("perl $exec_dir/get_gbk.pl  $basename  .");# create .gbk file
			last if (-e "$basename.gbk");
		}
		# we do not want bacterail genome to contaninate our virus db. If cannot download the .gbk file for 
		# virus identification, we dump this NC. 
		if (!-s "$basename.gbk"){
			print "  No $basename.gbk back from NCBI, skip\n";
			next;
		}
		if (`cat $basename.gbk` !~/\n\s*ORGANISM.*?(Viruses|Virus|phage).*?\nREFERENCE/is){
				#means the curr NC is not virus or phage.
			  print "    $basename is bacterial, skip\n";
				next;
		}
		$count++;
		print "     $v_name/$file\n";
		my $data = `cat $file`;
		$data =~s/\[.*?\]/\[$v_name\]/gs;
		$file =~s/\.faa//;
		$data =~s/>gi/>$prefix\_$file-gi/gs;
		print OUT $data;
		push @curr_arr, $file;
	}
	chdir $curr_dir;

}
close IN;
print "In all.faa.tar.gz, we complement $count entries to old virus!\n";
#now handle the not intercepted NC in phage.html
@need_arr = minus(\@need_arr, \@curr_arr);
my @old_keys = keys %hash_old_virus;
@need_arr = minus(\@need_arr, \@old_keys);
print "The rest virus genome count = ". (scalar @need_arr)." from phage.html\n";
# handle extra genomes
foreach my $s (@need_arr){
		#print "   $s from need_arr not found in curr_arr, call genbank to get the gbk file!\n";
		foreach (1..4){
    	system("perl $exec_dir/get_gbk.pl  $s  $tmp_dir");# create .gbk file
			last if (-s "$tmp_dir/$s.gbk");
		}
    system("perl $exec_dir/gbk2faa.pl $s.gbk"); #create .faa file
		if (-s "$tmp_dir/$s.gbk"){
			print "  $s.gbk back from genbank\n";
		}else{
		  print "  No $s.gbk back from GENbank\n";
		}
    my $data = `cat $s.faa`;
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
	      $l =~s/>gi/>$substr\_$s-gi/;
	  	}
	   print OUT $l."\n";
		}
		unlink "$tmp_dir/$s.gbk", "$tmp_dir/$s.faa";
}
close OUT;
my $error_flag = 0;
if (-s "virus.db"){
	print "append virus.db* to upper dir\n";
	system("cat virus.db >> ../virus.db");
	#cleanup 
	#system("rm -rf *") if (-s "list");
}else{
	print "virus.db  has no content. Please check!\n";
	$error_flag = 1;
}

# handle vgenome.tbl file for scan.pl to use
push @curr_arr,@need_arr; 
@curr_arr = uniq(@curr_arr) if (scalar @curr_arr !=0);
make_vgenome_file(@curr_arr);

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

sub make_vgenome_file{
	my @array = @_;
	    #check if different with the old phage list file
			my $data = `cat ../vgenome.tbl`;
			my @arr2 = ();
			foreach my $l (split("\n", $data)){
				my @tmp = split("\t", $l);
				push @arr2, $tmp[0];
			}
			my @rest = minus(\@array, \@arr2);
			if (scalar @rest ==0){
			   print "XXXXX There is no extra new virus entry from NCBI\n\n";
				 print "exit make vgenome file\n";
				 return;
			}else{
			   print "XXXXX There are ".(scalar @rest)." new virus entries from NCBI.\n\n";
			}
			open (OUT, ">vgenome.tbl") or die "Cannot append vgenome.tbl";
			my $c =0;
			foreach my $NC (@rest){
			  $c++;
			  print "\tWorking on $NC, count $c\n";
			  system("perl $exec_dir/get_gbk.pl  $NC  $tmp_dir");# create .gbk file
			  system("perl $exec_dir/gbk2faa.pl $tmp_dir/$NC.gbk"); #create .faa filed
			  my $data = `cat $tmp_dir/$NC.gbk`;
			  my $legth = '';
			  my $sp = '';
			  if ($data =~/LOCUS\s+\S+\s+(\d+)\s+bp.*?\nSOURCE\s+(.*?)\n/s){
			    $length=$1;
			    $sp = $2;
			    $sp =~s/\(.*\)//;
			  }
			  my $str ='';
			  my $count=0;
		    open (IN, "$tmp_dir/$NC.faa") or die "Cannot open $NC.faa";
		    while (my $l = <IN>){
		     if ($l=~/>gi\|(\d+)\|/){
		       $str .= "\t$1";
		       $count++;
			   }
		    }
		    close IN;
		    if ($count !=0){
		      print OUT "$NC\t$sp\t$length\t$count$str\n";
					print "          add $NC to vgenome.tbl\n";
	     	}else{
					print "          NO add $NC to vgenome.tbl, there is no #gi\n";
				}
		    #unlink "$tmp_dir/$NC.gbk", "$tmp_dir/$NC.faa";
		  }
		close OUT;
		print "Append vgenome.tbl to ../vgenome.tbl\n";
		system("cat vgenome.tbl >> ../vgenome.tbl");
		print "XXXX vgenome.tbl regenerated!\n";
}

sub get_prefix{
	my $v_name= shift;
	$v_name =~s/_phage//;
	$v_name =~s/Phage_//;	
	if ($v_name =~/^([A-Za-z]+)_/){
		my $str = $1;
		$str_sub = substr($str, 0, 6);
		$v_name =~s/^([A-Za-z]+)_/$str_sub\_/;
	}else{
		#print "Wierd, $v_name\n";
	}
	$v_name = "PHAGE_$v_name";		
	return $v_name;

}
sub uniq {
		my @a = @_;
		my @tmp=();
		my %hash=();
		foreach (@a){
		  next if (! defined $_ || $_ =~/^\s*$/);
			$hash{$_}=1;
		}
		@tmp = keys %hash;
		#@tmp = keys %{{ map { $_ => 1 } @a }} ;
		return @tmp;
}
sub minus{
 my $arr1= shift;
 my $arr2=shift;
 my @arr1 = (); @arr1 = uniq(@$arr1) if (scalar @$arr1 !=0);
 my @arr2 = (); @arr2 = uniq(@$arr2) if (scalar @$arr2 !=0);
 my @tmp=();
 foreach my $e1 (@arr1){
   my $f = 0;
   foreach my $e2 (@arr2){
     if ($e1 eq $e2){
       $f=1;
       last;
     }
   }
   push  @tmp, $e1  if ($f==0);
 }
 return @tmp;
}

sub get_hash_old_virus_db{
	my %hash =();
	open(IN, "../virus.db") or die "Cannot open ../virus.db";
	while(<IN>){
		if ($_=~/^>.*?_(\w{2}_?\d+)-gi\|(\d+)/){
			$hash{$1} = 1;
		}
	}
	close IN;
	return %hash;
}
