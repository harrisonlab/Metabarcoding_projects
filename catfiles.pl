#!/usr/bin/perl -w -s

my $its1 = $ARGV[0] || die "Please provide ITS1" ;
my $its2 = $ARGV[1] || die "Please provide ITS2" ;
my $myfasta = $ARGV[2] || die "Please provide output fasta" ;

my %hash_ITS1 = seqbuilder($its1);
my %hash_ITS2 = seqbuilder($its2);

open(MYFASTA, '>', $myfasta) or die "Could not open file '$myfasta' $!";

while( my( $key, $val1 ) = each( %hash_ITS1 ) ) {
	my( $val2 ) = $hash_ITS2{$key};
	if ($val1 ne "" && $val2 ne "")  {
		print MYFASTA "$key\n$val1$val2\n";
	}
}	

close(MYFASTA);

sub seqbuilder {
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
