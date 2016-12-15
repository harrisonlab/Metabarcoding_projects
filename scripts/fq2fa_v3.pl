#!/usr/bin/perl -s -w

# get filenames
my $usage = "$0 <sequence_file> <outfile> <sample ID> <trim_left> <trim_right> \n";
die $usage unless $ARGV[0];

#my @sequence_file = do { local $/; <STDIN> };

#my $sequence_file = $ARGV[0] || die "Please provide sequence file" ;
#my $outfile2 = $ARGV[0] || die "Please provide outfile name" ;
my $id = $ARGV[0] || die "Please provide sample ID" ;
my $ltrim=0;

if (defined $ARGV[1]) {
	$ltrim=$ARGV[1];	
}


#open (FILE, "<$sequence_file") || die "File $sequence_file doesn't exist!!";

# open output for appending
#open (OUTFILE2, ">>$outfile2") or die "Failed to open file '$outfile2' for writing\n";

# parse field and trim 
my $count = 2;
my $fid="";


while(<STDIN>) {
	if($_=~/\@/) {
		my $fas = $count/2;
		print ">$id.$fas;$_";
	}else {
		my $id_line=substr($_,$ltrim);
		print "$id_line";
	}
}

#while (my $id_line = <FILE>) {
#    if ($count%4==2) {
#    	$fid=$id_line;
#    }
#    $count++;
#    next if ($count%4!=0);
#    my $fas = $count/4;
#    print OUTFILE2 ">$id.$fas;$fid";
#    $id_line=substr($id_line,$ltrim,length($id_line)-($ltrim+$rtrim+1));
#    print OUTFILE2 "$id_line\n";
#}

#close FILE;
#close OUTFILE2;
