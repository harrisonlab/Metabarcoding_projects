#!/usr/bin/perl -w -s

##########################################################################################################
# Input is stdin for list of FASTQ read headers and command line arg for FASTQ file                      #
# Output is FASTA of reads matching stdin        		                             		 #
##########################################################################################################

# slurp stdin
my $headers = do { local $/; <STDIN> };


my ($inpFile) = $ARGV[0] || die "Please provide fastq input" ; # this is actually fasta now 
unless(open(INFILE, $inpFile) ) {
	print("Cannot open input file \"$inpFile\"\n\n");
     	exit;
}
my @fastq = <INFILE>;
chomp(@fastq);
close INFILE;

my $nskipper = 0;
my $counter = 4;

foreach(@fastq) {
	if ($counter%2==0) {
		my @new = split(";",$_,);
		if(index($headers,$new[1]) !=-1) {
			$nskipper = 2;
		}
	}
	if($nskipper) {
		print"$_\n";
		$nskipper--;
	}
	$counter++;
}
