#!/usr/bin/perl -s -w
##############################################################################
# makes a biom with taxa metadata from usearch output biom able and taxa.txt #
##############################################################################

my ($taxaFile) = $ARGV[0] || die "Please provide taxa file" ; 
unless(open(INFILE, $taxaFile) ) {
	print("Cannot open taxa file \"$taxaFile\"\n\n");
     	exit;
}
my @taxa = <INFILE>;
chomp(@taxa);
close INFILE;
my %seq=();
foreach(@taxa) {
 $_ =~ s/"//g;
 @t = split(/,/);
 $seq{$t[0]}="{\"taxonomy\": [\"k__$t[1]\", \"p__$t[3]\", \"c__$t[5]\", \"o__$t[7]\", \"f__$t[9]\", \"g__$t[11]\", \"s__$t[13]\"]}";
}  

my ($inpFile) = $ARGV[0] || die "Please provide biom table" ;
unless(open(INFILE, $inpFile) ) {
	print("Cannot open biom table \"$inpFile\"\n\n");
     	exit;
}
my @biom = <INFILE>;
chomp(@biom);
close INFILE;

foreach(@biom) {
  if(($OTU)=($_=~/.*(OTU\d*)",/)){
    $_=~s/null/$seq{$1}/;
  }
  print"$_\n"	;
}
