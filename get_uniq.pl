#! /usr/bin/perl -w -s

my %seq=();

my $fasta="";
my $lastcount=0;
while (<>) {
	if ($_=~/\>/) {

		if ($fasta ne "") {
			if (exists $seq{$fasta}) {
				$seq{$fasta}++;
			} else {
				$seq{$fasta}=1;
			}
		}

		$fasta="";		
	} else{
		chomp;
		$fasta.=$_;
	}
}

if ($fasta ne "") {
	if (exists $seq{$fasta}) {
		$seq{$fasta}++;
	} else {
		$seq{$fasta}=1;
	}
}


my $counter=1;
foreach my $key (keys %seq) {
	print">uniq.$counter;size=$seq{$key};\n$key\n";
	$counter++;
}
