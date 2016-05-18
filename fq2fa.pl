#!/usr/bin/perl

use strict;
use warnings ;

my $usage = "$0 <sequence_file> <outfile2>\n";
die $usage unless $ARGV[0];

my $sequence_file = $ARGV[0] || die "Please provide sequence file" ;
my $outfile2 = $ARGV[1] || die "Please provide outfile2 name" ;

open (FILE, "<$sequence_file") || die "File $sequence_file doesn't exist!!";


open (OUTFILE2, ">>$outfile2") or die "Failed to open file '$outfile2' for writing\n";

my $count = 2;
while (my $id_line = <FILE>) {
    $count++;
    next if ($count%4!=0);
    my $fas = $count/4;
    print OUTFILE2 ">$ARGV[2].$fas\n";
    print OUTFILE2 "$id_line";
}


close FILE;
close OUTFILE2;
