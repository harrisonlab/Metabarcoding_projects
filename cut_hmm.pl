#!/usr/bin/perl -w -s

my $hmm_file = $ARGV[0] || die "Please provide hmm file" ;
my $outdir = $ARGV[1] || die "Please provide output directory" ;
my $taxa = $ARGV[2] || die "Please provide taxa" ;
unless(-e $outdir or mkdir $outdir) {die "Unable to create $outdir";}
open (HMM, "<$hmm_file") || die "HMM $hmm_file doesn't exist!!";
my $outfile = "";
my $small_hmm ="";
while (my $line = <HMM>) {
	$small_hmm.=$line;
	chomp($line);
	if ($line =~/$taxa/) {
		$outfile=substr($line,index($line," ")+2);
	}

	if ($line=~/\/\//) {
		if ($outfile =~/$taxa/) {
			open OUTPUT, ">", "$outdir/$outfile.hmm" or die "Failed to open file '$outfile.hmm' for writing\n";
			print OUTPUT $small_hmm;
			close(OUTPUT);
		}
		$small_hmm="";
		$outfile="";
	}
	
}
