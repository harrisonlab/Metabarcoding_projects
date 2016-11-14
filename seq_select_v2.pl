#!/usr/bin/perl -s

##########################################################################################################
# Input is stdin for list of FASTQ read headers and command line arg for FASTQ file                      #
# Output is FASTA of reads matching stdin        		                             		 #
##########################################################################################################

# slurp stdin
#my $headers = do { local $/; <STDIN> };

my %seq=();

#my $fasta="";
#my $lastcount=0;
while ( <STDIN> ) {
	chomp;
	$seq{$_}="";
}

my ($inpFile) = $ARGV[0] || die "Please provide fasta input" ; 
unless(open(INFILE, $inpFile) ) {
	print("Cannot open input file \"$inpFile\"\n\n");
     	exit;
}
my @fastq = <INFILE>;
chomp(@fastq);
close INFILE;

my $nskipper = 0;
my $counter = 2;

foreach(@fastq) {
	if ($counter%2==0) {
		my @new = split(";",$_);
		my $test = $new[1];
		if (exists($seq{$test})) {
			$nskipper = 2;
		}
	}
	if($nskipper) {
		print"$_\n";
		$nskipper--;
	}
	$counter++;
}
