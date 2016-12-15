#!/usr/bin/perl -w -s

#save command line arguments to variables
my $its1 = $ARGV[0] || die "Please provide ITS1" ;
my $its2 = $ARGV[1] || die "Please provide ITS2" ;
my $myfasta1 = $ARGV[2] || die "Please provide output fasta file" ;
my $myfasta2 = $ARGV[3] || die "Please provide output fasta file" ;
my $nojoin  = 1;

#copy ITS1 and ITS2 to hashes
my %hash_ITS1 = seqbuilder($its1);
my %hash_ITS2 = seqbuilder($its2);

#open output files
open(myfasta1, '>', "$myfasta1.fa") or die "Could not open file '$myfasta1' $!";
open(myfasta2, '>', "$myfasta2.fa") or die "Could not open file '$myfasta1' $!";


# Keep count of reads for new fasta headers
my $count = 1;

my @U_ITS1 = ();


#find and print unique ITS2 
my @U_ITS2 = ();
foreach (keys %hash_ITS2) {
    if (exists $hash_ITS1{$_}) {
	print myfasta1 ">$myfasta1.$count\n$hash_ITS1{$_}\n";
    	print myfasta2  ">$myfasta2.$count\n$hash_ITS2{$_}\n";
    	$count++;
    }
}
	
close(myfasta1);
close(myfasta2);

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