#!/usr/bin/perl -w

# call glimmer parallel
use File::Basename;

my $fna_file = $ARGV[0];
my $log_file = "$fna_file.log";
my $cgi_dir = "/home/prion/phage/cgi";
my $dir = dirname($fna_file);
my $fna_file_basename = basename($fna_file);
my $offset=10000;
my $length_piece = 1000000; # the length of DNA sequence on each file .
my $data='';
open (IN, $fna_file) or die "Cannot open $fna_file";
my $header ='';
while (<IN>){
	if ($_=~/>/){
		$header = $_;
	}else{
		chomp($_);
		$data.=$_;
	}
}
close IN;
my $count = 0;
while(1){
	$count++;
	my $start;
	my $substr='';
	my $flag = 0;
	my $end_offset = 0;
	if ($count % 2==0){
		$start = $count/2;
		if (length($data) <= $start*$length_piece+ $offset){
			$substr= substr($data, $start*$length_piece-$offset);
			$end_offset = length($data)-$start*$length_piece;
			$flag = 1;
		}else{
			$substr= substr($data, $start*$length_piece-$offset, 2*$offset);
			$end_offset = $offset;
		}
		open(OUT, ">$fna_file\_$count") or die "Cannot write $fna_file\_$count";
		print OUT ">$count\t".($start*$length_piece-$offset+1)."\t".($start*$length_piece+$end_offset)."\n";
		print OUT $substr;
		close OUT;
		if ($flag ==1){
			last;
		}
	}else{
		$start = ($count-1)/2;
		if (length($data) <= $start*$length_piece+ $length_piece){
			$substr = substr($data, $start*$length_piece);
			$end_offset = length($data)-$start*$length_piece;
			$flag =1;
		}else{
			$substr = substr($data, $start*$length_piece, $length_piece);
			$end_offset = $length_piece;
		}
		open(OUT, ">$fna_file\_$count") or die "Cannot write $fna_file\_$count";
		print OUT ">$count\t".($start*$length_piece+1)."\t".($start*$length_piece+$end_offset)."\n";
		print OUT $substr;
		close OUT;
		if ($flag ==1){
			last;
		}
	}
}
system("echo 'piece $count' >> $log_file");

# Generate string specifying which nodes to run on.
my $last_node = 8; # max cluster node number to use
my $node_set = '';
for (my $i = 1; $i <= $last_node; $i++) { # restrict to particular nodes
	$node_set .= 'botha-w' . $i;
	if ($i != $last_node) {
		$node_set .= '|';
	}
}

# call glimmer pararllel
chdir $dir;
system("echo change to $dir >> $log_file");
system("echo qsub -t 1-$count   -q one.q -l h=\"$node_set\" -sync yes  $cgi_dir/single_glimmer.pl  $fna_file_basename >> $log_file");
system("qsub -t 1-$count -q one.q -l h=\"$node_set\" -sync yes  $cgi_dir/single_glimmer.pl  $fna_file_basename")==0 or system("echo $! >> $log_file");
open(OUT, ">$fna_file.predict") or die "Cannot write $fna_file.predict";
print OUT $header;
my $orf_count =0;
for (my $i=1; $i<=$count; $i++){
	if ($i % 2 ==0) {
		open(IN, "$fna_file_basename\_$i.predict") or die "Cannot open $fna_file_basename\_$i.predict";
		my $start_num;
		while (<IN>){
			if ($_=~/^>/){
				my @arr=split("\t",$_);
				$start_num = $arr[1];
				next;	
			}else{
				my @arr=split(" ", $_);
				$arr[1] += $start_num-1;
				$arr[2] += $start_num-1;
				if ($arr[2] >=$start_num-1+$offset && $arr[1] <=$start_num-1+$offset ) {
					$orf_count++; 
					my $orf_str='';
					if ($orf_count >0 && $orf_count<10) {
						$orf_str='0000'.$orf_count;
					}elsif($orf_count >=10 && $orf_count<100) {
						$orf_str='000'.$orf_count;
					}elsif($orf_count >=100 && $orf_count<1000) {
						$orf_str='00'.$orf_count;
					}elsif($orf_count >=1000 && $orf_count<10000) {
						$orf_str='0'.$orf_count;
					}else{
						$orf_str=$orf_count;
					}
					if (($arr[3]=~/\+/ && $arr[1] < $arr[2]) or ($arr[3]=~/\-/ && $arr[2] < $arr[1])) {
						
						my $line = sprintf("%-10s %-10s %-10s %-5s %-20s\n", "orf$orf_str", $arr[1], $arr[2], $arr[3], $arr[4]);
						#system("echo '$line' >>$log_file");
						print OUT $line;
					}
				}
				
			}
		}
		close IN;
	}else{
		open(IN, "$fna_file_basename\_$i.predict") or die "Cannot open $fna_file_basename\_$i.predict";
		my $start_num;
		while (<IN>){
			if ($_=~/^>/){
				my @arr=split("\t",$_);
				$start_num = $arr[1];
				next;	
			}else{
				my @arr=split(" ", $_);
				$arr[1] += $start_num-1;
				$arr[2] += $start_num-1;
				if (($arr[3]=~/\+/ && $arr[1] < $arr[2]) or ($arr[3]=~/\-/ && $arr[2] < $arr[1])) {
					$orf_count++; 
					my $orf_str='';
					if ($orf_count >0 && $orf_count<10) {
						$orf_str='0000'.$orf_count;
					}elsif($orf_count >=10 && $orf_count<100) {
						$orf_str='000'.$orf_count;
					}elsif($orf_count >=100 && $orf_count<1000) {
						$orf_str='00'.$orf_count;
					}elsif($orf_count >=1000 && $orf_count<10000) {
						$orf_str='0'.$orf_count;
					}else{
						$orf_str=$orf_count;
					}
					my $line = sprintf("%-10s %-10s %-10s %-5s %-20s\n", "orf$orf_str", $arr[1], $arr[2], $arr[3], $arr[4]);
					
					print OUT $line;
				}
				
			}
					
		}
		close IN;
	}
}
close OUT;
exit;


	
