#! /usr/bin/perl
use warnings;
use strict;

############################################################
#
# Demultiplex paired fastq with a set of regex indices
# Unidentified sequences are written to both output files
#
# Usage: demulti.pl forward_read reverse_read i1 i2
#
###########################################################

open (my $R1,'<', $ARGV[0]) or die "Could not open R1: $ARGV[1]";
open (my $R2,'<', $ARGV[1]) or die "Could not open R2: $ARGV[2]";

open (my $R1_bac,'>', "$ARGV[0].bacterial.fastq") or die "Could not open $ARGV[0].bacterial.fastq for writing";
open (my $R1_fun,'>', "$ARGV[0].fungal.fastq") or die "Could not open $ARGV[0].fungal.fastq for writing";
open (my $R2_bac,'>', "$ARGV[1].bacterial.fastq") or die "Could not open $ARGV[1].bacterial.fastq for writing";
open (my $R2_fun,'>', "$ARGV[1].fungal.fastq") or die "Could not open $ARGV[1].fungal.fastq for writing";

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
		if ($l1=~/$ARGV[2]/) {
			$output = 1;	
		} elsif ($l1=~/$ARGV[3]/) {
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
			print $R1_bac "$R1_l";
			print $R2_bac "$R2_l";
			print $R1_fun "$R1_l";
			print $R2_fun "$R2_l";
			$R1_l="";
			$R2_l="";
		}
	}
	$counter++;
} until eof $R1;
