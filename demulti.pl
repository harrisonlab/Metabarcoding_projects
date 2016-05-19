#! /usr/bin/perl -w -s

############################################################
#
# Demultiplex paired fastq with a set of regex indices
# Unidentified sequences are written to both output files
#
# Usage: demulti.pl forward_read reverse_read i1 i2 i3 i4 mismatch
#
###########################################################

open (my $R1_bac,'>', "$ARGV[0].bacterial.fastq") or die "Could not open $ARGV[0].bacterial.fastq for writing";
open (my $R1_fun,'>', "$ARGV[0].fungal.fastq") or die "Could not open $ARGV[0].fungal.fastq for writing";
open (my $R2_bac,'>', "$ARGV[1].bacterial.fastq") or die "Could not open $ARGV[1].bacterial.fastq for writing";
open (my $R2_fun,'>', "$ARGV[1].fungal.fastq") or die "Could not open $ARGV[1].fungal.fastq for writing";
open (my $R1_ambig,'>', "$ARGV[0].ambig.fastq") or die "Could not open $ARGV[1].ambig.fastq for writing";
open (my $R2_ambig,'>', "$ARGV[1].ambig.fastq") or die "Could not open $ARGV[1].ambig.fastq for writing";

open (my $R1,'<', shift) or die "Could not open R1: $ARGV[0]";
open (my $R2,'<', shift) or die "Could not open R2: $ARGV[1]";

my $ind_bac_f = shift;
my $ind_bac_r = shift;
my $ind_fun_f = shift;
my $ind_fun_r = shift;

my $allowmiss = shift;
$allowmiss ||= 1 unless defined $allowmiss;

my $output = 0;
my $R1_l = "";
my $R2_l = "";
my $counter=1;
do {
	my $l1 = <$R1>;
	my $l2 = <$R2>;

	$R1_l.=$l1;
	$R2_l.=$l2;

	if(($counter+2)%4==0) {

		my $countbf = count_match(substr($l1,0,length($ind_bac_f)),$ind_bac_f);
		my $countbr = count_match(substr($l2,0,length($ind_bac_r)),$ind_bac_r);
		my $countff = count_match(substr($l1,0,length($ind_fun_f)),$ind_fun_f);
		my $countfr = count_match(substr($l2,0,length($ind_fun_r)),$ind_fun_r);
		if ($countbf <= $allowmiss && $countbr <= $allowmiss) {
			$output = 1;
		} elsif ($countff <= $allowmiss && $countfr <= $allowmiss ) {
			$output = 2;

		} else {
			$output = 0;
		}
	}

	if($counter%4==0) {
		if($output==1) {
			print $R1_bac "$R1_l";
			print $R2_bac "$R2_l";
			$R1_l="";
			$R2_l="";
		} elsif ($output==2) {
			print $R1_fun "$R1_l";
			print $R2_fun "$R2_l";
			$R1_l="";
			$R2_l="";
		} else {
			#write to ambiguous sequence list
			print $R1_ambig "$R1_l";
			print $R2_ambig "$R2_l";
			$R1_l="";
			$R2_l="";
		}
	}

	$counter++;
} until eof $R1;


sub count_match {
	my $str1 = uc(shift);
	my $str2 = uc(shift);
	my $count = length($str1);
	for (my $i=0;$i<length($str1);$i++) {
		my $x = to_bin(substr($str1,$i,1));
		my $y = to_bin(substr($str2,$i,1));
		$count -- if $x&$y;
 	}
	return($count);
	
}

sub to_bin{
	my $t = @_[0];
	my $score = 0;
	$score ++ if $t=~/[ARWMDHVN]/;
	$score =$score+2 if $t=~/[CYSMBHVN]/;
	$score =$score+4 if $t=~/[GRSKBDVN]/;
	$score =$score+8 if $t=~/[TYWKBDHN]/;
	return($score);	
}
