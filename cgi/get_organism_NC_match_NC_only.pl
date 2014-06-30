#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  get_organism_NC_match.pl
#
#        USAGE:  ./get_organism_NC_match.pl  
#
#  DESCRIPTION:  This program will take the bacterial whole package's extract
#  		 list file as input and output the match file between organisms 
#  		 and NC numbers. it is possible to match one orgranism to many 
#  		 NC numbers. OUtput file name is ~/phage/DB/z_current_organism_NC_list.
# 		 another out file is the NC list only, ~/phage/DB/z_current_NC_list.
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Yongjie Liang(), 
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  06/05/2012 03:38:51 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
chdir "/home/prion/phage/DB/temp_bac_select";

my %hash=();
my @NC=();
open(IN, "zz_tar_list_file") or die "Cannot open zz_tar_list_file";
while (my $l=<IN>){
	my @a =split("/", $l);
	$a[1]=~s/_ui.+//;
	$a[2]=~s/\.faa\n//;

	push @NC, $a[2];
	unless(defined $hash{$a[1]}){
		$hash{$a[1]}='';	
	}
	$hash{$a[1]}.=$a[2].",";
}
close IN;
open(OUT, ">zz_current_organism_NC_list");
foreach my $k (sort {$a cmp $b} keys %hash){
	$hash{$k}=~s/,$//;
	print OUT "$k\t$hash{$k}\n";
}
close OUT;
open(OUT, ">zz_current_NC_list");
foreach my $k (sort {$a cmp $b} @NC){
	print OUT "$k\n";
}
close OUT;
system("mv zz_current_organism_NC_list ../z_current_organism_NC_list; mv zz_current_NC_list ../z_current_NC_list");
exit;


