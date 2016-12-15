#! /usr/bin/perl -w -s

my %seq=();

my $fasta="";
my $lastcount=0;
while (<>) {
	if ($_=~/\>/) {
		@l=split(/;size=/,$_);
		my $count = substr($l[1],0,-2);
		if ($fasta ne "") {
			if (exists $seq{$fasta}) {
				$seq{$fasta}+=$lastcount;
			} else {
				$seq{$fasta}=$lastcount;
			}
		}
		$lastcount=$count;
		$fasta="";		
	} else{
		chomp;
		$fasta.=$_;
	}
}

if ($fasta ne "") {
	if (exists $seq{$fasta}) {
		$seq{$fasta}+=$lastcount;
	} else {
		$seq{$fasta}=$lastcount;
	}
}


my $counter=1;
foreach my $key (keys %seq) {
	print">uniq.$counter;size=$seq{$key};\n$key\n";
	$counter++;
}