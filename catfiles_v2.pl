#!/usr/bin/perl -w -s

#save command line arguments to variables
my $its1 = $ARGV[0] || die "Please provide ITS1" ;
my $its2 = $ARGV[1] || die "Please provide ITS2" ;
my $myfasta = $ARGV[2] || die "Please provide output fasta file" ;
my $nojoin  = 1;

#copy ITS1 and ITS2 to hashes
my %hash_ITS1 = seqbuilder($its1);
my %hash_ITS2 = seqbuilder($its2);

#open output files
open(MYFASTA, '>', "$myfasta.fa") or die "Could not open file '$myfasta' $!";


#Print ITS1 
my @U_ITS1 = ();
foreach (keys %hash_ITS1) {
    print MYFASTA "$_\n$hash_ITS1{$_}\n";
}

#find and print unique ITS2 
my @U_ITS2 = ();
foreach (keys %hash_ITS2) {
    print MYFASTA "$_\n$hash_ITS2{$_}\n" unless exists $hash_ITS1{$_};
}
	
close(MYFASTA);

sub seqbuilder {
# Returns hash of input fasta file 
# Key is fasta header and value is sequence

	my ($inpFile) = @_;
    	unless(open(INFILE, $inpFile) ) {
		print("Cannot open input file \"$inpFile\"\n\n");
        	exit;
    	}
	my @file = <INFILE>;
	chomp(@file);
	close INFILE;
	my %sequences = ();
	my $label;
	foreach(@file) {
		chomp;
		if ($_=~m/^\>/){
			$label = $_;
			$sequences{$label}="";
		}
		else{
			$sequences{$label}.=uc($_);
		}
	}
	return %sequences;
}