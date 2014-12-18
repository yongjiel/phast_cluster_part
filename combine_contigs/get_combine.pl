#!/usr/bin/perl -w 
#===============================================================================
#
#         FILE:  get_combine.pl
#
#        USAGE:  ./get_combine.pl  
#
#  DESCRIPTION:  This program will get gbk contig files and combine them into a
#                single genome.
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Jack Liang
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  04/03/2013 12:16:04 PM
#     REVISION:  ---
#===============================================================================

use File::Path;
use File::Copy;

if (@ARGV != 1){
	print STDERR "Usage: perl get_combine.pl <basename, like: NZ_ALJA00000000> \n";
	exit(-1);
}

my $basename = $ARGV[0];
if (-s "../result_tmp/$basename/$basename.gbk"){
	print "$basename is in result/tmp already\n";
	exit(0);
}
my $tmp_dir = "tmp_$basename";
  rmtree("$tmp_dir");
	  mkdir "$tmp_dir";
system("perl get_gbk.pl $basename $tmp_dir")==0 or print "$!\n"; 
if (!-s "$tmp_dir/$basename.gbk"){
	print STDERR "Cannot get $basename.gbk from genbank\n";
	exit(-1);
}
my ($start, $end)=`cat $tmp_dir/$basename.gbk`=~/\nWGS\s+(.*?)-(.*?)\n/s;
if (!defined $start or !defined $end){
	($start, $end)=`cat $tmp_dir/$basename.gbk`=~/\nWGS_SCAFLD\s+(.*?)-(.*?)\n/s;
}
if (!defined $start){
	print STDERR "Undefine start in $tmp_dir/$basename.gbk\n";
	exit(-1);
}
if (!defined $end){
	print STDERR "Undefine end in $tmp_dir/$basename.gbk\n";
	exit(-1);
}

my $basename1=$start; $basename1=~s/(\d{6})$/000000/;
my $from =int($1);
my ($to) =$end; $to=int($1) if ($to=~/(\d{6})$/);
my $last_len=0;
my @array = ();
for (my $i =$from ; $i <= $to; $i++){
	my $s='';
	if(length($i) < length($to)){
		for(my $j=1; $j <= length($to)-length($i); $j++){
			$s.='0';
		}
		$s.=$i;
	}else{
		$s = $i;
	}
	my $len= length($s);
	my $single_file = $basename1; $single_file =~s/\d{$len}$/$s/;
	if (!-s "$tmp_dir/$single_file.gbk"){
		system("perl get_gbk.pl $single_file  $tmp_dir  ")==0 or print "$!\n";
	}
	if (!-s "$tmp_dir/$single_file.gbk"){
		print STDERR "Something Wrong. No $tmp_dir/$single_file.gbk\n";
		exit(-1);
	}else{
		push @array,  "$single_file.gbk";
	}

}
my $gbk_file= "$basename.gbk";
open (OUT, "> $gbk_file");
my $seq='';
foreach $f (@array){
	open (IN, "$tmp_dir/$f");
	my $start_flag = 0;
	my $seq_flag = 0;
	my $sequ='';
	while (my $line =<IN>){
		if ($line=~/^\s+gene\s+/){
					$start_flag = 1;
		}
		if ($line=~/^ORIGIN/){
					$start_flag = 0;
					$seq_flag = 1;
					next;
		}
		if ($line=~/^\/\//){
			$seq_flag = 0;
		}
		if ($start_flag ==1){
			$line=~s/(\d+)\.\.(.*?)(\d+)/($1+$last_len)."..$2".($3+$last_len)/ge;
			print OUT $line;
		}
		if ($seq_flag ==1){
			$line=~s/[\s\d\n]//g;
			$sequ.=$line;
			$seq.=$line;
		}
	}
	close IN;
	my $head = `head -n 1 $tmp_dir/$f`;
	my ($l)= $head =~/(\d+) bp/;
	if (!defined $l or $l eq ''){
		print STDERR "$l not defined \n";
		exit(-1);
	}
	if ($l != length($sequ)){
		print "$l != ".(length($sequ))."in $tmp_dir/$f\n";
	}
	$last_len +=$l;
	print $last_len."\n";
}
close OUT;

my ($def_part) =`cat $tmp_dir/$basename.gbk`=~/(\nDEFINITION\s+.*?\nVERSION\s+.*?\n)/s;
my $CDS_part = `cat $gbk_file`;
my $first_file = "$tmp_dir/$array[0]";
my $f_p = `cat $first_file`;
$f_p=~s/\nORIGIN.*//s;
$f_p=~s/\n\s+gene\s+.*//s;
my $v = '00000000';
$f_p=~s/LOCUS       ([^\d]+)\d+/LOCUS       $1$v/s;
$f_p=~s/\d+ bp/$last_len bp/s;
$f_p=~s/\n\s+source\s+1..\d+/\n     source          1..$last_len/s;
$f_p=~s/bases 1 to \d+/bases 1 to $last_len/gs;
$f_p=~s/\nDEFINITION\s+.*?\nVERSION\s+.*?\n/$def_part/s;
my @seq = $seq=~/\w{0,10}/g;
my $seq_str = '';
my $n = 1;
my $count =0;
foreach my $s (@seq){
	next if($s eq '');
	if ($count==0){
		$seq_str .= sprintf("%9s", $n);
	}
	$count++;
	if ($count %6==0 ){
		$seq_str.=" ".$s."\n";
		$n+=60;
		$seq_str .= sprintf("%9s", $n);
	}else{
		$seq_str .= " $s";
	}
}
open (OUT, ">$gbk_file");
print OUT $f_p."\n";
print OUT $CDS_part;
print OUT "ORIGIN\n";
print OUT $seq_str;
close OUT;

rmtree $tmp_dir;
mkdir "../result_tmp/$basename";
copy  $gbk_file, "../result_tmp/$basename/$gbk_file"; 
system("ssh -i /home/prion/.ssh/scp-key phast\@phast.wishartlab.com 'echo \"-a $basename\" >> ~/project/tmp/queue.txt' ");
exit;

