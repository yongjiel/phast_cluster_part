#!/usr/bin/perl -w

use Bio::Perl;
use Bio::DB::GenBank;

#Usage: perl get_gbk.pl  <accession_number or gi number>  <tmp_dir>
# <flag> is -g or -a ; -g means using gi number, and -a using accession number

my $dir = $ARGV[1];
my $flag = '';
if ($ARGV[0]=~/^\d+$/){
	$flag = '-g';
}else{
	$flag = '-a';
}

$gb = Bio::DB::GenBank->new();
if ($flag eq '-a'){
	$seq= $gb->get_Seq_by_acc($ARGV[0]);
}
elsif ($flag eq '-g'){
	$seq= $gb->get_Seq_by_gi($ARGV[0]);
}

$out = Bio::SeqIO->new(-file => ">$dir/$ARGV[0].gbk", -format => 'Genbank');
$out->write_seq($seq);

exit;

