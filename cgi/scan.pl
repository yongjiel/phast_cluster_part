#!/usr/bin/perl
# April 9, 2011 Rah. This script process preannoated gbk inputs.
# This will make PHAST faster to handle preprocessed input cases because blast only runs on
# identified (by protein names) regions (instead of the entire genome)
# Note: additional work will be needed for remove transposons
# example command
#perl phast.pl -g NC_002662/NC_002662.gbk -n NC_002662/NC_002662.fna -a NC_002662/NC_002662.faa -t NC_002662/tRNAscan.out -m NC_002662/tmRNA_aragorn.out -p NC_002662/NC_002662.ptt -b NC_002662/ncbi.out

use strict;
use lib "/apps/phast/source";
use Gene;
use Dictionary;
use Prophage;
use PhageTable;
use Record;
use DNA;
# server specific variable
my $logdir = "/apps/phast/log/";
my $dbdir = "/apps/phast/database/";
my $sub_program_dir ="/apps/phast/sub_programs";
# genbank file alone is sufficient for all operations but current script uses all files
my $fna = "";
my $trn = "";
my $tmr = "";
my $faa = "";
my $bla = ""; # blast result (ncbi.out), please use "null" if not available
my $ptt = "";
my $gbk = ""; # genbank file
my $log = $logdir."phast.log"; # log file

# steps:
# 1. identify all phast regions by keywords (words related to phage):
# 2. run local blast on identified regions against the phage database
# or use returned results of whole genome blast if given

my $seedRadius = 3;		# half size of the mandatory blast region (each side of the key words)
my $extentionWindow = 10;	# extention window size
my $minRate = 0.3;			# min phage hit rates to put the new CDS to blast queue

# 3. assess the blasted regions to find prophages
# 4. identified integrases and (pre-identified)tRNAs are used to predict att sites
# 5. bacterial protein information is copied from gbk file to identified prophages
# 6. completeness is assigned for each prophage using phage genome information
# 7. predicted prophage are printed in phage finder's format

# blast parameters:
my $tophit = 5;			# number of hits returned
my $evalue = 0.001;		# cut-off e value of hit
my $processors = 1;		# number of processors to use
my $bprogram = "blastall"; 		# blast program
my $db = "$sub_program_dir/phage_finder/DB/prophage_virus.db";	# blast database
my $usetop = 250;		# max top hits to use

# DBSCAN parameters:
my $eps = 3000;			# distance parameter (see DBSCAN paper)
my $minpts = 4;			# minimal CDS to form a phage

# premature prophages joining parameters
my $maxIntDistance = 10000;	# this is the max distance an "noise" integrase can be joined to the nearby prophage (otherwise the integrase will be made into an individual prophage)
my $phageGenomeTable = $dbdir."vgenome.tbl";
my $maxphage = 150000;		# maximal phage in bps after join

# attachment site search
# some parameters are coded in DNA.pm
#my $minlen = 12; # see DNA.pm
#my $maxlen = 25;

# database and tables
#my $virustable = $dbdir."phage_db_anno.tbl";	# points to a table of virus annotation
my $virustable = `grep '>' $db`;
# time stamp
my $ProgramStartTime = time;

# debug var
my $debug_p_afterDBSCAN = 0;
my $debug_p_afterJoin = 0;
my $debug_p_afterblast = 0;
my $filter = 1;

# parse arguments
if ($#ARGV == -1){
	&help();
}
for my $i (0 .. $#ARGV){
	if ($ARGV[$i] eq '-h'){
		&help();
	}
	elsif ($ARGV[$i] eq '-g'){
		$gbk = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-n'){
		$fna = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-a'){
		$faa = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-t'){
		$trn = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-b'){
		$bla = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-m'){
		$tmr = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-p'){
		$ptt = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-use'){
		$usetop = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-Deps'){
		$eps = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-Dpts'){
		$minpts = $ARGV[++$i];
	}
	elsif ($ARGV[$i] eq '-dDBSCAN'){
		$debug_p_afterDBSCAN = 1;
	}
	elsif ($ARGV[$i] eq '-djoin'){
		$debug_p_afterJoin = 1;
	}
	elsif ($ARGV[$i] eq '-dblast'){
		$debug_p_afterblast = 1;
	}
	elsif ($ARGV[$i] eq '-rmfilter'){
		$filter = 0;
	}
	else{
#		die "unknown flag '$ARGV[$i] at position' $i";
	}
}

# dynamic data
my $sequence = new DNA($fna);
my $genomehead = ''; 	# head of sorted Gene list
my $genometail = ''; 	# tail of sorted Gene list
my $prophagehead = '';		# head of sorted prophage region list
my %localgihash = ();	# hash table keyed by Gene local gi number (exclude RNA)
my %starthash = ();		# hash table keyed by Gene start index (include RNA)


# local functions
sub help {
	print "PHAST v1.0 Dr. Whishart's lab (http://phast.wishartlab.com)\n";
	print "Please contact Rah Zhou (email: youz\@ualberta.ca)\n";
	print "flags:";
	print "  -h                 print this help page\n";
	print "  -g [file]          Genbak file (.gbk)\n";
	print "  -n [file]          DNA sequence file (.fna)\n";
	print "  -a [file]          amino acid (protein) sequence file (.faa)\n";
	print "  -t [file]          tRNAscan-SE output file\n";
	print "  -m [file]          Aragorn output file\n";
	print "  -b [file]          NCBI blast result (ncbi.out, tuple format -8)\n";
	print "  -p [file]          genome information file (.ptt)\n";
	print "  -Btop		        specify the top number of blast results for analysis\n";
	print "  -Beva              specify blast evalue\n";
	print "  -Deps              specify eps parameter of DBSCAN clustering algorithm\n";
	print "  -Dpts              specify minimal number cluster size of DBSCAN clustering algorithm\n";
	print "  -B                 do not use blast result\n";
	print "In order to run either Genbank file or .ptt file has to be given.\nIf Genbank file is given annotations will be used for identifying phage genes in the genome.\n";
	exit;
}
# insert a Gene to this sorted list. This could be expensive because the list can be long
sub insertGene {
	my $gene = shift;
	my $cursor;
	if ($genomehead eq ''){
		$genomehead = $gene;
		$genometail = $gene;
		return 0;
	}
	elsif ($gene->getStart() <= $genomehead->getStart()){
		$gene->setNext($genomehead);
		$genomehead->setPrevious($gene);
		$genomehead = $gene;
	}
	elsif ($gene->getStart() >= $genometail->getStart()){
		$gene->setPrevious($genometail);
		$genometail->setNext($gene);
		$genometail = $gene;
	}
	elsif (abs($gene->getStart() - $genomehead->getStart()) < abs($gene->getStart() - $genometail->getStart())){
		for ($cursor = $genomehead; $cursor->getStart() < $gene->getStart(); $cursor = $cursor->next()){}
		$gene->setNext($cursor);
		$gene->setPrevious($cursor->previous());
		$cursor->previous()->setNext($gene);
		$cursor->setPrevious($gene);
	}
	else{
		for ($cursor = $genometail; $cursor->getStart() > $gene->getStart(); $cursor = $cursor->previous()){}
		$gene->setPrevious($cursor);
		$gene->setNext($cursor->next());
		$cursor->next()->setPrevious($gene);
		$cursor->setNext($gene);
	}
}

sub findCoverlapRNA {
	my ($start, $end) = @_;
	my $cursor;
	for ($cursor = $genomehead; $cursor ne ''; $cursor = $cursor->next()){
		if ($cursor->isRNA()){
			if ($cursor->getEnd() >= $start  and $cursor->getStart() <= $end){
				return $cursor;
			}
		}
		elsif ($cursor->getStart() > $end){
			last;
		}
	}
	return '';
}

sub findGeneByLocalGI {
	my $targetGI = shift;
	return $localgihash{$targetGI};
}

sub findGeneByStart {
	my $targetStart = shift;
	return $starthash{$targetStart};
}

# genbank file alone is sufficient for all operations but current script uses all files
# future work will be to write this function to use gbk file only to reduce complexity of the code
sub buildGenome{
	# parse .ptt file (assume proteins are sorted by index)
	open (R, $ptt) or die "cannot open ptt file $ptt, ptt is mandatory";
	while (my $line = <R>){
	# Location        Strand  Length  PID     Gene    Synonym Code    COG     Product
	# 19..1677        +       552     38232643        dnaA    DIP0001 -    -   chromosomal replication initiation protein
		if ($line =~m/^(\d+)\.\.(\d+)\t(.*?)\t.*?\t(.*?)\t.*?\t(.*?)\t.*?\t.*?\t(.*)/){
			my $gene = new Gene($4, $1, $2, $3, $6, 'p', $5); # local gi, start, end, strand, product and type
			&insertGene($gene);
			$localgihash{$4} = $gene;
			$starthash{$1} = $gene;
		}
		elsif ($line =~m/^\d+\.\.\d+/){
			print STDERR "DIE: error parsing ptt file at line: $line\n";
		}
	}
	close (R);
	# parser .gbk file to insert tRNA (if exists)
	if (-e $gbk){
		open (R, $gbk) or die "cannot read genbank file $gbk (hardware)";
		while (my $line = <R>){
			my $strand = '';
			my $start = '';
			my $end = '';
			my $product = '';
			if ($line=~m/^\s+tRNA\s+<*>*(\d+)\.+<*>*(\d+)/){
				$strand = '+';
				$start = $1;
				$end = $2;
			}
			elsif ($line=~m/^\s+tRNA\s+complement\(<*>*(\d+)\.+<*>*(\d+)/){
				$strand = '-';
				$start = $1;
				$end = $2;
			}
			if ($start ne ''){
				while ($line = <R>){ # looking for tRNA product
                                        if ($line=~m/^\s\s\s\s\s\w/){last;} # move on ifhit the next feature (assume no consective tRNAs)
                                        if ($line=~m/^\s+\/product="(.+)"/){
                                                $product = $1;
                                                last;
                                        }
                                }
				#print "$start, $end, $strand, $product, 'r'\n";
				my $gene = new Gene('', $start, $end, $strand, $product, 'r'); # local gi, start, end, strand, product and type
				&insertGene($gene);
				$starthash{$start} = $gene;
			}
		}
		close (R);
	}
	#print "start parse tRNA file\n";
	if (-e $trn){
		open (R, $trn) or die "cannot read $trn which should be the tRNA file (hardware)";
		while (my $line = <R>){
			# Name                        tRNA #	Begin  	End    	Type	Codon	Begin	End	Score
			# gi|16271976|ref|NC_000907.1| 	1	47152  	47241  	Ser	TGA	0	0	69.75
			if ($line =~m/^gi/){
				my @tokens = split (/\t/, $line);
				if ($#tokens < 8){
					print STDERR "error: tRNA file line: $line\n";
				}
				my ($strand, $start, $end);
				if ($tokens[2] < $tokens[3]){
					$strand = '+';
					$start = $tokens[2];
					$end = $tokens[3];
				}
				else{
					$strand = '-';
					$start = $tokens[3];
					$end = $tokens[2];
				}
				# local gi, start, end, strand, product and type
				my $tRNAgene = new Gene('', $start, $end, $strand, "tRNA-$tokens[4]", 'r');
				my $existRNA =  &findGeneByStart($start);
				if ($existRNA eq ''){
					$existRNA = &findCoverlapRNA($start, $end);
				}
				if ($existRNA eq ''){
					&insertGene($tRNAgene);
				}
			}
		}
		close (R);
	}
	#print "start parse tmRNA\n";
	if (-e $tmr){
		open (R, $tmr) or die "cannot read $trn which should be the tmRNA file (hardware)";
		while (my $line = <R>){
			# seems can either be "Location c[1357564,1357929]" or "Location [1357564,1357929]"
			my ($strand, $start, $end);
			if ($line =~m/^Location\s*\[(\d+),(\d+)\]/){
				$strand = '+';
				$start = $1;
				$end = $2;
			}
			elsif ($line =~m/^Location\s*c\[(\d+),(\d+)\]/){
				$strand = '-';
				$start = $1;
				$end = $2;
			}
			else{
				next;
			}
			# local gi, start, end, strand, product and type
			my $tRNAgene = new Gene('', $start, $end, $strand, "tmRNA", 'r');
			my $existRNA =  &findGeneByStart($start);
			if ($existRNA eq ''){
				$existRNA = &findCoverlapRNA($start, $end);
			}
			if ($existRNA eq ''){
				&insertGene($tRNAgene);
			}
		}
		close (R);
	}
	#print "start parse (ncbi.out) blast file\n";
	if (-e $bla){
		open (R, $bla) or die "cannot read $bla, which is the viral database blast result (hardware)";
		my %blast = (); # store results
		my @all = <R>;
		close (R);
		foreach(@all){
			my @tokens = split(/\t/, $_);
			if ($tokens[0]=~m/^gi\|(\d+)/){
				my $querygi = $1;
				#PHAGE_Prochl_P-SSM2-gi|61806059|ref|YP_214419.1|
				if ($tokens[1]=~m/gi\|(\d+)\|ref\|(.+)\|/){
					my $gi = $1; 
					my $ref = $2;
					my $sp = '';
					if ($tokens[1]=~m/^(.*)-gi\|\d+\|ref\|.+\|/){
						$sp = $1;	
					}
					my $record = new Record($gi, $sp, "", $ref, $tokens[10]); # Gi, species and definition
					if (defined $blast{$querygi}){
						my $array =  $blast{$querygi};
						if ($#$array+1 < $usetop){
							$$array[$#$array+1] = $record;
						}
					}
					else{
						my @array = ();
						push (@array, $record);
						$blast{$querygi} = \@array;
					}
				}
				else{
					print STDERR "WRONG FORMAT in blast out: $_";
				}
			}
		}
		while (my ($k, $v) = each %blast){
			my $target = &findGeneByLocalGI($k);
			if ($target ne ''){
				$target->setBLASTresult($v);
			}
			else{
				print STDERR  "error: blast returned unexpected results for local GI $k\n";
				next;
			}
		}
	}
}

# print the genome list (protein and rRNA) for debug
sub printGenome {
	my ($hits) = @_;
	for (my $cursor = $genomehead; $cursor ne ''; $cursor = $cursor->next()){
		$cursor->printInfo($hits);
	} 
}

# print the prophage list for debug
sub printProphage {
	my ($hits) = @_;
	for (my $cursor = $prophagehead; $cursor ne ''; $cursor = $cursor->next()){
		$cursor->printInfo($hits);
	}
}

# system call to blast against phage database
sub blast {
	# pick up CDS from queue
	my %gi = ();
	for (my $cursor = $genomehead; $cursor ne ''; $cursor = $cursor->next()){
		if ($cursor->getBLASTstatus() == 1){
			$gi{$cursor->getLocalGI()} = $cursor;
			$cursor->setBLASTstatus(2); # mark as blasted
		}
	}
	# make the .faa file for blast
	open (R, $faa) or die "cannot open .faa file";
	my @all = <R>;
	close (R);
	my $tfile = "b".time;
	open (W, ">$tfile.bi") or die "cannot write temperatory file";
	for (my $i = 0; $i <= $#all; $i++){
		if ($all[$i]=~m/^>gi\|(\d+)\|/){
			if (defined $gi{$1}){
				print W $all[$i++];
				if ($i > $#all){
					last;
				}
				elsif ($all[$i]=~m/^\w+/){
					print W $all[$i++];
				}
				$i--;
			}
		}
	}
	close (W);
	# call blast
	system ("$bprogram -p blastp -d $db -m 8 -e $evalue -i $tfile.bi -o $tfile.bo -v 4 -b $tophit -a $processors -F F");
	# parse returned file
	open (R, "$tfile.bo") or die "cannot open blast output $tfile.bo";
	@all=<R>;
	close (R);
	system ("rm $tfile.b*");
	my %blast = (); # store results
	foreach(@all){
		my @tokens = split(/\t/, $_);
		if ($tokens[0]=~m/^gi\|(\d+)/){
			my $querygi = $1;
			#print "GI $1";
			if ($tokens[1]=~m/^(.*[^-])-?gi\|(\d+)/){
				#print " GI $2  NAME $1\n";
				if (defined $blast{$querygi}){
					my $array =  $blast{$querygi};
					$$array[$#$array+1] = "$2|$1";
				}
				else{
					my @array = ();
					push (@array, "$2|$1");
					$blast{$querygi} = \@array;
				}
			}
			else{
				print STDERR  "WRONG FORMAT in blast out!$tokens[1]\n";
			}
		}
	}
	while (my ($k, $v) = each %blast){
		if (defined $gi{$k}){
			$gi{$k}->setBLASTresult($v);
			#print "set $k\n";foreach(@$v){print "$_\n";} 
		}
		else{
			print STDERR "error: blast returned unexpected results for GI $k\n";
		}
	}
	print "One BLAST run finished\n";
}

# sub to annotate the blast result (because blast resut has no annotation)
sub addAnnotation {
	#open (R, $virustable) or die "cannot read virus annotation table at: $virustable";
	my %hash = ();
	foreach my $line (split("\n", $virustable)){
	  #while (my $line = <R>){
		chomp ($line);
		#if ($line=~m/^(\d+)\t(.+)/){
		#>PROPHAGE_Xantho_306-gi|21243357|ref|NP_642939.1| phage-related integrase [Xanthomonas axonopodis pv. citri str. 306]
		if ($line=~/gi\|(\d+)\|ref\|.*?\|\s*(.*)/){
			$hash{$1} = $2;
		}
	}
	#close (R);
	for (my $cursor = $genomehead; $cursor ne ''; $cursor = $cursor->next()){
		my $hits = $cursor->getBLASTresult();
		if ($hits ne '' and $hits ne "keyword"){
			foreach (@$hits){
				if (defined $hash{$_->getGI()}){
					$_->setDefinition($hash{$_->getGI()});
				}
				else{
					$_->setDefinition("ERROR IN DATABASE: CANNOT LINK THIS PROTEIN");
					print STDERR "ERROR: virus database has no record for GI: ", $_->getGI(), "\n";
				}
			}
		}
	}
}

# prophage list
sub insertProphage {
	my $prophage = shift;
	if ($prophagehead eq ''){
		$prophagehead = $prophage;
	}
	elsif ($prophage->getHead()->getStart() < $prophagehead->getHead()->getStart()){
		$prophage->setNext($prophagehead);
		$prophagehead->setPrevious($prophage);
		$prophagehead = $prophage;
	}
	else{	
		for (my $cursor = $prophagehead; $cursor ne ''; $cursor = $cursor->next()){
			if ($cursor->next() eq ''){
				$prophage->setPrevious($cursor);
				$cursor->setNext($prophage);
				last;
			}
			elsif ($prophage->getHead()->getStart() < $cursor->next()->getHead()->getStart()){
				$prophage->setNext($cursor->next());
				$prophage->setPrevious($cursor);
				$cursor->next()->setPrevious($prophage);
				$cursor->setNext($prophage);
				last;
			}
		}
	}
}

sub numberOfProphage {
	my $c = 0;
	for (my $cursor = $prophagehead; $cursor ne ''; $cursor = $cursor->next()){
		$c++;
	}
	return $c;
}
# DBSCAN algorithms
#DBSCAN(D, eps, MinPts)
#   C = 0
#   for each unvisited point P in dataset D
#      mark P as visited
#      N = getNeighbors (P, eps)
#      if sizeof(N) < MinPts
#         mark P as NOISE
#      else
#         C = next cluster
#         expandCluster(P, N, C, eps, MinPts)
#          
#expandCluster(P, N, C, eps, MinPts)
#   add P to cluster C
#   for each point P' in N 
#      if P' is not visited
#         mark P' as visited
#         N' = getNeighbors(P', eps)
#         if sizeof(N') >= MinPts
#            N = N joined with N'
#      if P' is not yet member of any cluster
#         add P' to cluster C
#end
sub DBSCAN_main {
	my @noise = ();
	for (my $unvisite = $genomehead; $unvisite ne ''; $unvisite = $unvisite->next()){
		if (!$unvisite->isRNA() and $unvisite->getBLASTresult() ne '' and $unvisite->isVisited() == 0){
			$unvisite->setVisited(1);
			my $neighbors = &DBSCAN_neighbors($unvisite);
			if ($#$neighbors+1 < $minpts){
				push (@noise, $unvisite);
			}
			else{
				my $prophage = new Prophage();
				&DBSCAN_expand($unvisite, $neighbors, $prophage);
				&insertProphage($prophage);
			}
		}
	}
	
	# savage any integrase from noise
	# there're usually integrase in noise array, we try to put them to closest
	# prophages
	foreach (@noise){
		if (&isPhageGene($_) == 0 and $_->isIntegrase() == 1){
#			$_->printInfo();
		# find the previous and next phage neighbor $_
			my $pre = '';
			my $nex = '';
			for ($nex = $prophagehead; $nex ne ''; $nex = $nex->next()){
				if ($nex->getHead()->getStart()>$_->getStart()){
					$pre = $nex->previous();
					last;
				}
			}
			if ($prophagehead eq ''){
				# no prophage found yet
			}
			elsif ($pre eq '' and $nex ne ''){
				if ($nex->getHead()->getStart - $_->getEnd() < $maxIntDistance){
					$nex->insert($_);
#					print "add to the first phage\n";
					next;
				}
			}
			elsif ($pre eq '' and $nex eq ''){
				for ($pre = $prophagehead; $pre->next() ne ''; $pre = $pre->next()){}
				if ($_->getStart() - $pre->getTail()->getEnd() < $maxIntDistance){
					$pre->insert($_);
#					print "add to the last phage\n";
					next;
				}
			}
			else{
				if ($nex->getHead()->getStart - $_->getEnd() < $maxIntDistance 
				and $_->getStart() - $pre->getTail()->getEnd() > $maxIntDistance){
					$nex->insert($_);
#					print "add to the right phage because left is out of range\n";
					next;
				}
				elsif ($nex->getHead()->getStart - $_->getEnd() > $maxIntDistance 
				and $_->getStart() - $pre->getTail()->getEnd() < $maxIntDistance){
					$pre->insert($_);
#					print "add to the left phage because right is out of range\n";
					next;
				}
				elsif ($nex->getHead()->getStart - $_->getEnd() < $maxIntDistance 
				and $_->getStart() - $pre->getTail()->getEnd() < $maxIntDistance){
					if ($pre->hasIntegrase() == 1 and $nex->hasIntegrase() == 0){
						$nex->insert($_);
#						print "add to the right phage because left has int but not right\n";
						next;
					}
					elsif ($pre->hasIntegrase() == 0 and $nex->hasIntegrase() == 1){
						$pre->insert($_);
#						print "add to the left phage because right has int but not left\n";
						next;
					}
					else{
						if ($_->getStart() - $pre->getTail()->getStart() < $nex->getHead()->getStart() - $_->getEnd()){
							$pre->insert($_);
#							print "add to the left phage because its closer\n";
							next;
						}
						else{
							$nex->insert($_);
#							print "add to the right phage because its closer\n";
							next;
						}
					}
				}
			}
#	April 12, 2011 by Rah: commented out since they seem to introduce false positives			
#			my $neighbors = &DBSCAN_neighbors($_);
#			my $prophage = new Prophage();
#			$prophage->insert($_);
#			for my $i (0 .. $#$neighbors){
#				$prophage->insert($_);
#			}
#			print "inserted as an independent phage\n";
#			&insertProphage($prophage);
		}
	}
}

sub DBSCAN_neighbors {
	my $self = shift;
	my @neighbors = ();
	for (my $cursor = $self->previous(); $cursor ne '' and $self->getStart() - $cursor->getEnd() < $eps; $cursor = $cursor->previous()){
		if ($cursor->getBLASTresult() ne ''){
			push (@neighbors, $cursor);
		}
	}
	for (my $cursor = $self->next(); $cursor ne '' and $cursor->getStart() - $self->getEnd() < $eps; $cursor = $cursor->next()){
		if ($cursor->getBLASTresult() ne ''){
			push (@neighbors, $cursor);
		}
	}
	return \@neighbors;
}

sub DBSCAN_expand {
	my ($p, $N, $c) = @_; # meanings of names see psudocode
	$c->insert($p);
	for (my $i = 0; $i <= $#$N; $i++){
		if ($$N[$i]->isVisited == 0){
			$$N[$i]->setVisited(1);
			my $M = &DBSCAN_neighbors($$N[$i]);
			if ($#$M+1 >= $minpts){
				foreach (@$M){
					push (@$N, $_);
				}
			}
		}
		if (&isPhageGene($$N[$i]) == 0){
			$c->insert($$N[$i]);
		}
	}
}

sub structNode {
	my %phage = ();
	my $self = {
		_phage => \%phage,
		_score => 0,
		_penalty => 0,
		_next => ''
	};
	return $self;
}
sub joinProphage {
	# join prophages if they are found to have the same phage genome
	my $table = new PhageTable($phageGenomeTable, $filter);
	# check all possible combination of prophages region for better results
#	printProphage(1);
	my @toptable = ();
	for (my $px = $prophagehead; $px ne ''; $px = $px->next()){
		my $penalty = 0;
		my $headnode = '';
		for (my $py = $px; $py ne '' and $py->getTail()->getEnd() - $px->getHead()->getStart() < $maxphage; $py = $py->next()){
			my $node = &structNode();
			if ($py == $px){
				$penalty = 0;
				$headnode = $node;
				push (@toptable, $headnode);
				$table->clear();
#				print "HEAD at ($py)",$py->getHead()->getStart(), "-", $py->getTail()->getEnd()," (P:$penalty)\n";
			}
			else{
				$headnode->{_next} = $node;
				$headnode = $node;
				for (my $cur = $px; $cur != $py; $cur = $cur->next()){
					$node->{_phage}->{$cur} = 1;
				}
				$penalty += $py->getHead()->getStart() - $py->previous()->getTail()->getEnd();
#				print "    extend at ",$py->getHead()->getStart(), "-", $py->getTail()->getEnd()," (P:$penalty)\n";
			}
			$table->assessprophage($py);
			$node->{_phage}->{$py} = $py;
			$node->{_score} = $table->evaluate(1000000, $penalty, $maxphage);
			$node->{_penalty} = $penalty;
#			if ($px == $py){
#				print "START: ", $px->getHead()->getStart(),"\n";
#				print $table->printParts();
#				print ">>>", $py->getHead()->getStart()," score = $node->{_score}\n";
#				#exit;
#			}
#			print "        score = $node->{_score}\n";
		}
	}
#	exit;
	# pick-up best prophages
	my @best = ();
	my %picked = ();
#	print "pick up best phage\n";
	while (1){
#		print "######### NEW RUN ########\n";
		my $p = 0;
		my $bestnode = &structNode();	# a node with 0 score
		foreach (@toptable){
#			print "Table Head: $_\n";
			for (my $node = $_; $node ne ''; $node = $node->{_next}){
				my $removed = 0;
#				print "  Node $node  $node->{_phage}\n";
				while (my ($phage, $v) = each %{$node->{_phage}}){
#					print "    contains prophage $phage ";
					if (defined $picked{$phage}){
#						print "ALREDY REMOVED\n";
						$removed = 1;
			#			last; # must allow while loop complete or the hash iterator won't reset 
					}
#					print "may be removed\n";
				}
				if ($removed == 1){
					last;
				}
				elsif ($node->{_score} > $bestnode->{_score}){
#					print "  score $node->{_score} > $bestnode->{_score} updated\n";
					$bestnode = $node
				}
				$p++;
			}
		}
#		print "best score ($p) = ", $bestnode->{_score},"\n";
		if ($bestnode->{_score} >= 1000000 - 10000){
			push (@best, $bestnode);
			while (my ($phage, $v) = each %{$bestnode->{_phage}}){
#				print "add choosen phage $phage ", $bestnode->{_score},"\n";
				if (defined $picked{$phage}){print STDERR "$phage: reached unexpected situation\n";}
				$picked{$phage} = $bestnode;
			}
		}
		else{
			last;
		}
	}
	
#	print "rebuild prophage list\n";
	# rebuild the prophage list
	my $newhead = '';
	foreach (@best){
#		print "BEST $_\n";
		my $newphage = new Prophage();
		while (my ($phage, $v) = each %{$_->{_phage}}){
#			print "  add phage $phage to  $newphage\n";
			my $ptr;
			# put the class reference to a hash and take it back. Perl won't recognize it, has to do
			# this trick to find it by value.
			for ($ptr = $prophagehead; $ptr ne ''; $ptr = $ptr->next()){
				if ($ptr eq $phage){
					last;
				}
			}
			
			for (my $gene = $ptr->getHead(); $gene ne ''; $gene = $gene->next()){
				$newphage->insert($gene);
				if ($gene == $ptr->getTail()){
					last;
				}
			}
		}
		if ($newhead eq ''){
			$newhead = $newphage;
		}
		else{
			$newphage->setNext($newhead);
			$newhead = $newphage;
		}
	}
#	print "out of while\n";
	$prophagehead = '';
	my $cursor = $newhead;
	while ($cursor ne ''){
		my $toAdd = $cursor;
		$cursor = $cursor->next();
		$toAdd->setNext('');
		$toAdd->setPrevious('');
#		print "insert $toAdd\n";
		&insertProphage($toAdd);
	}
	
	# April 24, 2011. Add completeness to Prophage
	for ($cursor = $prophagehead; $cursor ne ''; $cursor = $cursor->next()){
		$table->clear();
		$table->assessprophage($cursor);
		$cursor->setCompleteness($table->completeness());
	}
}

# attachment site
sub searchAttachmentSite {
#	my $t = $prophagehead;
#	for (my $i = 0; $i < 0; $i++){$t = $t->next();}
#	$sequence->findAttachmentSite($t); 
	for (my $cur = $prophagehead; $cur ne ''; $cur = $cur->next()){
#		print "Prophage: $cur\n";
		$sequence->findAttachmentSite($cur); # this also sets att in prophage
#		print "END $cur\n";
		$cur->calculateGC($sequence);	# add GC content information (after p5, p3 determinated)
	}
	
}
# subs to handle prophage list

sub isPhageGene {
	my $gene = shift;
	for (my $cursor = $prophagehead; $cursor ne ''; $cursor = $cursor->next()){
		if ($cursor->isInclude($gene)){
			return 1;
		}
	}
	return 0;
}

# print to STDOUT output in Phage Finder's format
sub printoutput {
	open (W, ">>$log") or die "cannot write to log file $log";
	# header:
	my $asml = "NC_000000";
	if ($faa =~m/(NC_\d+)/){
		$asml = $1;
	}
	print W "START:$asml:$ProgramStartTime\n";
	my $pindex = 1;
	for (my $p = $prophagehead; $p ne ''; $p = $p->next(), $pindex++){
		print $sequence->getDefninition();
		print " [asmbl_id: $asml].";
		print $sequence->getSize(), ", ";
		print sprintf("gc\%: %5.2f\%\n", 100.0*$sequence->getGlobalGC()/$sequence->getSize());
		print "Medium degenerate region $pindex is from ", $p->get5end(),
		" to ", $p->get3end()," and is ", $p->getSize()," bp in size, ",
		sprintf("gc\%: %5.2f\%", $p->calculateGC($sequence)), ".", $p->getCompleteness(),"\n";
		$p->printPhagefinderResult();
		print W "PHAGE:$pindex:", $p->get5end(),":", $p->get3end(),"\n";
	}
	print W "END:$asml:$ProgramStartTime:", time - $ProgramStartTime,"\n";
	close W;
	print sprintf ("There are %d regions between 10 and 18 Kb and summing 26486 bp of sequence ( 1.45\% of the genome)\n", $pindex-1);
}
##########  main program ###########

# reusable variables
my ($cursor, $cursor1, $cursor2, $i);

&buildGenome(); # load input file to make the genome list
#&printGenome();



if (! -e $bla){ # no blast (ncbi.out) given, do local blast use Dr. Wishart's instruction
	die "no ncbi.out given '$bla'";
	# step 1. identify all phast regions by keywords:
	for ($cursor = $genomehead; $cursor ne ''; $cursor = $cursor->next()){
		if ($cursor->isRNA() == 0 and &generalRelated($cursor->getProduct())){ # sub in Dictionary.pm
			$cursor->setRelationToPhage(1);
			print $cursor->getStart(),"  '",$cursor->getProduct(),"'\n";
		}
	}
	# step 2. run local blast on identified regions against the phage database
	# there's a practical problem how to choose close CDS for blast. This problem is addressed
	# using following logic:
	# 1: n CDS of each side and each the name identified CDS are mandatory blast targets
	# 2: a window of m CDS first moves from left to right, if number of blast hits exceed x percent then add
	#    the new CDS from the right to blast queue (if it has not been blasted)
	# 3: similar window moves from right to left.
	# 4: process terminates if when no new CDS is added to the queue

	for ($cursor = $genomehead; $cursor ne ''; $cursor = $cursor->next()){
		if ($cursor->isRelatedToPhage()){
			$i = 0;
			for ($cursor1 = $cursor; $i < $seedRadius and $cursor1 ne ''; $cursor1 = $cursor1->previous(), $i++){
				if (!$cursor1->isRNA()){
					$cursor1->setBLASTstatus(1); # mark to run blast
				}
			}
			$i = 0;
			for ($cursor1 = $cursor; $i < $seedRadius and $cursor1 ne ''; $cursor1 = $cursor1->next(), $i++){
				if (!$cursor1->isRNA()){
					$cursor1->setBLASTstatus(1); # mark to run blast
				}
			}
		}
	}
	&blast();

	my $newblast = 0; # note: number may double count
	do{
		$newblast = 0;
		my $hit = 0; # number of phage related CDS in the window
		# 5->3 window
		for ($i = 0, $cursor = $genomehead; $i < $extentionWindow and $cursor ne ''; $cursor = $cursor->next()){
			if (!$cursor->isRNA()){
				if ($cursor->getBLASTresult() ne ''){
					$hit++;
				}
				$i++;
			}
		}
		# cursor is the new CDS added and cursor1 is the last CDS removed from the window
		for ($cursor = $cursor->next(), $cursor1 = $genomehead;
			$cursor ne '';
			$cursor = $cursor->next(), $cursor1 = $cursor1->next()){
			# this assume cursor has been set to correct position
			if (!$cursor->isRNA()){
				if ($cursor->getBLASTresult() ne ''){
					$hit++;
				}
				while ($cursor1->isRNA()){
					$cursor1 = $cursor1->next();
				}
				if ($cursor->getBLASTresult() ne ''){
					$hit--;
				}
			}
			if ($hit / $extentionWindow >= $minRate and $cursor->getBLASTstatus() != 2){ # 2 means it has been blasted
				$cursor->setBLASTstatus(1);
				$newblast++;
			}
		}

		# 3->5 window
		$hit = 0;
		for ($i = 0, $cursor = $genometail; $i < $extentionWindow and $cursor ne ''; $cursor = $cursor->previous()){
			if (!$cursor->isRNA()){
				if ($cursor->getBLASTresult() ne ''){
					$hit++;
				}
				$i++;
			}
		}
		# cursor is the new CDS added and cursor1 is the last CDS removed from the window
		for ($cursor = $cursor->previous(), $cursor1 = $genometail;
			$cursor ne '';
			$cursor = $cursor->previous(), $cursor1 = $cursor1->previous()){
			# this assume cursor has been set to correct position
			if (!$cursor->isRNA()){
				if ($cursor->getBLASTresult() ne ''){
					$hit++;
				}
				while ($cursor1->isRNA()){
					$cursor1 = $cursor1->previous();
				}
				if ($cursor->getBLASTresult() ne ''){
					$hit--;
				}
			}
			if ($hit / $extentionWindow >= $minRate and $cursor->getBLASTstatus() != 2){ # 2 means it has been blasted
				$cursor->setBLASTstatus(1);
				$newblast++;
			}
		}
		if ($newblast > 0){
			&blast();
		}
	}while ($newblast > 0);
}
else{
	# blast has been pre-performed for all CDS, only need to identify phage by keywords in genbank file
	# identify prophages by names in genbank file
	# step 1 and 2 (done by read blast result). identify all phast regions by keywords:
	for ($cursor = $genomehead; $cursor ne ''; $cursor = $cursor->next()){
		if ($cursor->isRNA() == 0 and $cursor->getBLASTresult() eq '' and &strictRelated($cursor->getProduct()) == 1){ # sub in Dictionary.pm
			$cursor->setRelationToPhage(1);
			if ($cursor->getBLASTresult() eq ''){
				$cursor->setBLASTresult("keyword");
			}
#			print $cursor->getStart(),"  '",$cursor->getProduct(),"'\n";
		}
	}
}

&addAnnotation();

if ($debug_p_afterblast == 1){
	&printGenome(1);
	exit;
}
# step 3. assess the blasted regions to find prophages (by DBSCAN)
&DBSCAN_main();
if ($debug_p_afterDBSCAN == 1){
	&printProphage(1);
	exit;
}

# step 3.5. join prophage regions if they seem to have considerable proteins from the same phage genome
&joinProphage();

if ($debug_p_afterJoin == 1){
	&printProphage(1);
	exit;
}
# step 4. identified integrases and (pre-identified)tRNAs are used to predict att sites
&searchAttachmentSite();

# finally, format output in Phage Finder's format
&printoutput();

exit;
