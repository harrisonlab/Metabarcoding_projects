#! /usr/bin/perl 
#use warnings;
use strict;


open (my $R1_bac,'>', "$ARGV[0].bacterial.fastq") or die "Could not open $ARGV[0].bacterial.fastq for writing";
open (my $R1_fun,'>', "$ARGV[0].fungal.fastq") or die "Could not open $ARGV[0].fungal.fastq for writing";
open (my $R2_bac,'>', "$ARGV[1].bacterial.fastq") or die "Could not open $ARGV[1].bacterial.fastq for writing";
open (my $R2_fun,'>', "$ARGV[1].fungal.fastq") or die "Could not open $ARGV[1].fungal.fastq for writing";
open (my $R1_ambig,'>', "$ARGV[0].ambig.fastq") or die "Could not open $ARGV[1].ambig.fastq for writing";
open (my $R2_ambig,'>', "$ARGV[1].ambig.fastq") or die "Could not open $ARGV[1].ambig.fastq for writing";

open (my $r1,'<', shift) or die "Could not open R1: $ARGV[0]";
open (my $r2,'<', shift) or die "Could not open R2: $ARGV[1]";
open (my $m1,'<', shift) or die "Could not open R1: $ARGV[2]";
open (my $m2,'<', shift) or die "Could not open R2: $ARGV[3]";



my $lbf = length($ind_bac_f);
my $lbr = length($ind_bac_r);
my $lff = length($ind_fun_f);
my $lfr = length($ind_fun_r);
