#!/usr/bin/perl

use strict;
use warnings ;

# get filenames
my $usage = "$0 <sequence_file> <outfile> <sample ID> <trim_left> <trim_right> \n";
die $usage unless $ARGV[0];

my $sequence_file = $ARGV[0] || die "Please provide sequence file" ;
my $outfile2 = $ARGV[1] || die "Please provide outfile name" ;
my $id = $ARGV[2] || die "Please provide sample ID" ;
my $ltrim=0;
my $rtrim=0;

if (defined $ARGV[3]) {
	$ltrim=$ARGV[3];	
}

if (defined $ARGV[4]) {
	$rtrim=$ARGV[4];	
}

open (FILE, "<$sequence_file") || die "File $sequence_file doesn't exist!!";

# open output for appending
open (OUTFILE2, ">>$outfile2") or die "Failed to open file '$outfile2' for writing\n";

# parse field and trim 
my $count = 2;
while (my $id_line = <FILE>) {
    $count++;
    next if ($count%4!=0);
    my $fas = $count/4;
    print OUTFILE2 ">$id.$fas\n";
    $id_line=substr($id_line,$ltrim,length($id_line)-($ltrim+$rtrim+1));
    print OUTFILE2 "$id_line\n";
}


close FILE;
close OUTFILE2;
