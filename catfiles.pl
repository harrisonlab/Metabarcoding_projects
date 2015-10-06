#!/usr/bin/perl -w -s

#save command line arguments to variables
my $its1 = $ARGV[0] || die "Please provide ITS1" ;
my $its2 = $ARGV[1] || die "Please provide ITS2" ;
my $myfasta = $ARGV[2] || die "Please provide output fasta file" ;

#copy ITS1 and ITS2 to hashes
my %hash_ITS1 = seqbuilder($its1);
my %hash_ITS2 = seqbuilder($its2);

#open output files
open(MYFASTA, '>', "$myfasta.fa") or die "Could not open file '$myfasta' $!";
open(MYFASTAR1, '>', "$myfasta.r1.fa") or die "Could not open file '$myfasta.r1' $!";
open(MYFASTAR2, '>', "$myfasta.r2.fa") or die "Could not open file '$myfasta.r2' $!";

#find common fasta headers
foreach (keys %hash_ITS1) {
	if (exists $hash_ITS2{$_}) {
		my $val1=$hash_ITS1{$_};
		my $val2=$hash_ITS2{$_};
		print MYFASTA "$_\n$val1$val2\n";
	}
}

#find unique ITS1 fasta headers
my @U_ITS1 = ();
foreach (keys %hash_ITS1) {
    print MYFASTAR1 "$_\n$hash_ITS1{$_}\n" unless exists $hash_ITS2{$_};
}

#find unique ITS2 fasta headers
my @U_ITS2 = ();
foreach (keys %hash_ITS2) {
    print MYFASTAR2 "$_\n$hash_ITS2{$_}\n" unless exists $hash_ITS1{$_};
}

# This bit of code can be used if ITS1 and ITS2 have same fasta headers as alternative to finding common/unique headers (will save time for very long fasta files)
#while( my( $key, $val1 ) = each( %hash_ITS1 ) ) {
#	my( $val2 ) = $hash_ITS2{$key};
#	if ($val1 ne "" && $val2 ne "")  {
#		print MYFASTA "$key\n$val1$val2\n";
#	}
#}	

close(MYFASTA);
close(MYFASTAR1);
close(MYFASTAR2);

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

