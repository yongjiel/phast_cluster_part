#!/usr/bin/perl


my $input=$ARGV[0];
my $organismutput=$input;
$organismutput =~s/\.gbk//;
$organismutput.=".faa";
gbk2faa($input,$organismutput);

exit;
	
sub gbk2faa{
	my ($input,$output)=@_;
	my $organism=''; #organism
	my $gi='NONE'; #db_xref:"GI
	my $protein_id='NONE'; #protein_id
	my $product='NONE'; #product, protein name
	my $seq='NONE'; #amino acid sequence
	my $location ='NONE'; # coplement(114..567) pr 114..567

	my $gbk='';
	my @draft=();

	# Get the GenBank data into an array from a file
	$gbk = `cat $input`;
	$gbk =~s/.*?\n(FEATURES\s+.*?)\nORIGIN.*/$1/s;	
	if($gbk=~/\/organism="(.*?)"/s) {
		$organism=$1;
		$organism=~s/\n//g;
		$organism=~s/\[.*?\]//g;
		$organism=~s/\s\s+/ /g;
		$organism=~s/\s+$//;
        }
	$gbk =~s/^FEATURES\s+.*?(\s+gene\s+.*)/$1/s;
	
	@lines = split("\n", $gbk);
	my $line='';
	my $flag=0;
	
	foreach(@lines){
		if ($_=~/^\s*$/){
			next;
		}
		if ($_=~/^\s\s\s+gene\s\s\s+/){
			$flag=0;
			next;
		}
		if ($_=~/^\s\s\s+CDS\s\s\s+/){
			if ($line ne ''){
				push @draft, $line;
			}
			$line = $_."\n";
			$flag=1;
			next;
		}
		if ($flag==1){
			$line .= "$_\n";
		}
		
	}
	if ($line ne ''){
		push @draft, $line;
	}
	
	open(OUT,">$output");
	foreach my $rec (@draft){
		$gi='NONE';
		$protein_id='NONE';
		$product='NONE';
		$seq='NONE';
		$location ='NONE';
		
		if($rec=~/product="(.*?)"/s) {
			$product=$1;
			$product=~s/\n\s+/ /gs;	
        	}
		
		if ($rec=~/\s+CDS\s+(.*?)\//s ) {
			$location = $1;
			$location =~s/\n//g;
			$location =~s/\s//g;
			$location =~s/>//g;
			$location =~s/<//g;
			$location =~s/join\((\d+)\.\..*?(\d+)\)/$1\.\.$2/;
		}
		if($rec=~/protein_id="(.*?)"/s) {
			$protein_id=$1;
        	}
		if($rec=~/db_xref="GI:(.*?)"/s) {
			$gi=$1;
        	}
		if($rec=~/translation="(.*?)"/s) {
			$seq = $1;
			$seq =~s/[\s\n]//sg;
		}
		if ($seq eq "NONE"){
			next;
		}
		if ($gi eq 'NONE' or $protein_id eq 'NONE' or $product eq 'NONE' or $seq eq 'NONE' or $location eq 'NONE'){
			print "gi=$gi, protein_id=$protein_id, product=$product, seq=$seq, location=$location\n";
		}
		print OUT ">gi|$gi|ref|$protein_id| $product, $location [$organism]\n$seq\n";
	}
			
	close OUT;
	return 1;
}


