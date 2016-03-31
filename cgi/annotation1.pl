#!/usr/bin/perl -w


# Since there is no annotation in phage finder's output if we submit the raw nucleotide sequence 
# only into the server, now we complememt this part after phage finder finishes. In phage finder's 
# result, there is gi number there if there is blast hit on that region. We use that gi to get back
# the annotation from GENBank
# Usage : annotation.pl 

print "\nCall annotation.pl <NC number>  <dir-name>  <virus_db> <bac_db> <non_hit_blast_output_file> <scan_output_file> [-s|-a|-g]\n";
my $t1 = time();
my $NC=$ARGV[0];
my $num =$ARGV[1];
our $virus_db=$ARGV[2];
our $bac_db= $ARGV[3];
my $blast_output_file = $ARGV[4];
my $scan_output_file = $ARGV[5];
my $flag= $ARGV[$#ARGV];
my $blast_data='';

if (!-e $virus_db){
	print STDERR "There is no $virus_db, check\n";
	exit(-1);
}
if (!-e $bac_db ){
	print STDERR "There is no $bac_db , check\n";
	exit(-1);
}

my %hash_ptt_data=();
open(IN, "../$num.ptt") or die "Cannot open ../$num.ptt";
while(<IN>){
	if ($_=~/(\d+)\.\.\d+\t\+\t\S+\t(\S+)/ or $_=~/\d+\.\.(\d+)\t\-\t\S+\t(\S+)/){
		$hash_ptt_data{$1}= $2;
	}
}
close IN;
my $tt = time();
open(IN, "$virus_db") or die "Cannot open $virus_db";
my %hash_virus_database = ();
while(<IN>){
	if ($_=~/gi\|(\d+)\|\w{3}\|\S+\|\s*(.*?)\n/){
		$hash_virus_database{$1} = $2;
	}
}
close IN;
print "read in $virus_db, time ". (time()-$tt). " seconds\n";
my %hash_non_hit_blast_data=();

open(BAC_DB_IN, "$bac_db") or die "Cannot open $bac_db";
my %hash_faa_data=();
if ($NC eq 'NC_000000' or $flag eq '-s'){
	if (-s $blast_output_file){
		open(IN1, $blast_output_file) or die "Failed to open file: $!\n";
		while(<IN1>) { 
			chomp; 
			if ($_=~/^gi\|(\d+)\|.*?gi\|(\d+)\|.*?\s+(\S+)\s+[\d\.]+$/){
				$hash_non_hit_blast_data{$1}="$2 $3" if (!defined $hash_non_hit_blast_data{$1});
			}
		} 
		close IN1;
	}
}else{
	open(IN, "../$num.faa") or die "Cannot open ../$num.faa";
	#>gi|16273379|ref|NP_439626.1| hypothetical protein, complement(1559375..1559734) [Haemophilus influenzae Rd KW20]
	while(<IN>){
		if ($_=~/>gi\|(\d+)\|ref\|\S+\|(.*?)\n/s){
			my $gi=$1;
			$hash_faa_data{$gi} = $2;
			$hash_faa_data{$gi}=~s/complement\(<*>*\d+\.\.<*>*\d+\)//;
			$hash_faa_data{$gi}=~s/<*>*\d+\.\.<*>*\d+//;
			$hash_faa_data{$gi}=~s/\[.*?\]//;
			$hash_faa_data{$gi}=~s/,.*$//;
		}
	}
	close IN;
}

if (-e $scan_output_file){
	get_annotation($scan_output_file, \%hash_non_hit_blast_data, $num, \%hash_virus_database, \%hash_ptt_data, \*BAC_DB_IN, \%hash_faa_data);
	print "  Change $scan_output_file \n";
}
close BAC_DB_IN;
my $diff_time = time()- $t1;
print "Finish annotation.pl in $diff_time seconds\n\n";
exit(0);

# get annotation on specified file.
sub get_annotation{
	my ($scan_output_file ,$hash_non_hit_blast_data, $num, $hash_virus_database, $hash_ptt_data, $bac_db_fh, $hash_faa_hash)=@_;
	#my $title = `cat ../$num.fna`;
	my $title = '';
	open(TI, "<../$num.fna") or die "Cannot open file $!";
	while(<TI>){
		$title .= $_;
	}
   	close TI;
	$title =~s/>(gi.*?)\n.*/$1/s;
	$title =~s/[,\.]\s*$//;
	system("cp $scan_output_file $scan_output_file\_bk");
	my $phage_suffix='';
	$phage_suffix='phage' if ($scan_output_file=~/phmedio.txt/);
	open (IN, $scan_output_file);
	open (OUT, ">$scan_output_file.tmp") or die "Cannot write $scan_output_file.tmp";
	print "Writing $scan_output_file.tmp\n";
	my $desc;
	while(<IN>){
		if ($_=~/^\d+\s+\d+\s+gi\|(\d+)\|/){ # 1739269  01576              gi|9634162|ref|NP_037720.1|, TAG = PHAGE_coliph_HK97, E-VALUE = 1e-60
			print OUT $_;
			my $gi = $1;
			$desc = 'N/A';
			if (defined $hash_virus_database->{$gi}){
				$desc=$hash_virus_database->{$gi}; $desc=~s/\s*\[.*?\].*//;
				$desc=~s/<*>*\d+\.\.<*>*\d+//;
		#		print "hit virus db $gi\n";
			}
			my $next_line='';
			my $count =0;
			# sometimes it will look like this:
			#1758638  01611              gi|9635602|ref|NP_061499.1|, TAG = PHAGE_pseudo_D3, E-VALUE = 8e-07
			#gi|00000000|ref|NC_000000| Pseudomonas putida KT2440 chromosome, complete genome [asmbl_id: NC_000000], 6181863, gc%: 61.52%
			#Large prophage region 1 is from 1738082 to 1776397 and is 38316 bp in size, gc%: 61.10%.

			#END5     FEATNAME           ANNOTATION OF BEST HIT FROM PHAGE DB, [ANNO]TATION, OR [HMM] HIT                                                                       PAGE: 2
			#............................................................................................................................................................................
 			#                           [ANNO] -; pp_01611

			$next_line = <IN>; $desc =~s/\n//g;
			$next_line=~s/\[ANNO\].*?;\s*(\S+)/\[ANNO\] $desc; $1; $phage_suffix/; # [ANNO] -; pp_01713  ---> [ANNO] integrase, putative; pp_1532; phage
			#print $next_line."\n";
			#exit;
			print OUT $next_line;
		}elsif ($_=~/^(\d+)\s+\d+\s+\[ANNO\]\s*-.*?;\s*([\w_]+\d+)/i){ #1742858  01583              # [ANNO] -; pp_01713 
								#571689     00538              [ANNO] -  -; PP_00538
								#571689     00538              [ANNO] -  -; PP_00538
			my $e=''; my $gi='';
			my $end5=$1;
			
			my $gi_b = '';
			if (defined $hash_ptt_data->{$end5}){
			# 4579..5391      -       270     16271981        fdhD    HI0005  -       formate dehydrogenase accessory protein
				$gi_b = $hash_ptt_data->{$end5};
			}else{
				print "No $end5 in ptt_data\n";
			}	
			
			if ($NC eq 'NC_000000' or $flag eq '-s'){# raw DNA sequence
				if (defined $hash_non_hit_blast_data->{$gi_b}){
					($gi,$e) = $hash_non_hit_blast_data->{$gi_b}=~/(\d+) (\S+)/;
				}
			}else{ # have gbk file and faa file is from gbk .
				$gi = $gi_b;
				$e = 0.00;
			}
			my $line = $_;
			if ($gi ne ''){
				$desc = 'N/A';
				if ($NC eq 'NC_000000' or $flag eq '-s' ){
				#	if(defined $hash_bac_database->{$gi}){
				#		$desc=$hash_bac_database->{$gi};
				#	}
					$desc = get_decs_from_bac_db($bac_db_fh, $gi);
				}else{
					#>gi|16273379|ref|NP_439626.1| hypothetical protein, complement(1559375..1559734) [Haemophilus influenzae Rd KW20]
					if (defined $hash_faa_data->{$gi}){
						$desc = $hash_faa_data->{$gi};
					}
				}
				
				$desc =~s/\.$//; $desc =~s/\n//g;
				$line=~s/\[ANNO\]\s*-.*?;\s*([\w_]+\d+)/\[ANNO\] $desc; E-VALUE = $e; $1/i; #  # [ANNO] -; pp_01713 ---> [ANNO] protein; pp_1542; hypothetical
			}else{# gi is ''
				#967646     01049              [ANNO] -; BASYS_01049
				$line=~s/\[ANNO\]\s*-.*?;\s*([\w_]+\d+)/\[ANNO\] hypothetical; $1/i; #  # [ANNO] -; pp_01713 ---> [ANNO] hypothetical; pp_1542
			}
			#print "Changed\n";
			print OUT $line;
		}elsif ($_=~/^.*?\[.*?\].*gc\%/ && $_!~/gi/ ){ # [asmbl_id: NC_000000], , gc%: %
			$_= $title.$_;
			print OUT $_;
		}else{
			print OUT $_;
		}
	}
	close IN;
	close OUT;
	open(IN, "$scan_output_file.tmp") or die "Cannot open $scan_output_file.tmp";
	open(OUT, ">$scan_output_file") or die "Cannot write $scan_output_file";
	while (<IN>){
		print OUT $_;
	}
	close IN;
	close OUT;
	unlink "$scan_output_file.tmp";
}

sub get_decs_from_bac_db{
	my ($fh, $gi) = @_;
       my $t1 = time;
        my $rc = `grep '^$gi' $bac_db\.ind -A1`;
        my ($start_index, $end_index) = $rc =~ /^\d+\s+(\d+)\n\d+\s+(\d+)/s;
#        print "time = ".(time-$t1)."\n";
        $t1 = time;
 #       print "start_index, $start_index, end_index, $end_index\n";
        if (defined $start_index and defined $end_index){
     		my $len = $end_index - $start_index;
      		my $str;
      		seek($fh, $start_index, 0);
      		read($fh, $str, $len);
               #>gi|158333234|ref|YP_001514406.1| NUDIX hydrolase [Acaryochloris marina MBIC11017]
 		if ($str =~ m/>(gi\|\d+\|ref\|.*?\|)\s+(.*)/s){
			my $h = $1;
			my $desc = $2; $desc=~s/\n//;
	#		print "$desc. $h\n";
  #   print "time = ". (time - $t1) ."\n";
			return "$desc. $h";
                }else{
                        die "Cannot return desc for $gi!\n";
		}
        }else{
		die "Cannot return desc for $gi!\n";
	}
}


