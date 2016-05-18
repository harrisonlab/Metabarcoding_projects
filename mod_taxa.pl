#! /usr/bin/perl -s

use Switch;

my $d ="unknown,100";
my $p ="unknown,100";
my $c ="unknown,100";
my $o ="unknown,100";
my $f ="unknown,100";
my $g ="unknown,100";
my $s ="unknown,100";

while(<>) {

	my @l = split(/\t/,$_);
	test_it(($l[6],$l[7],$l[5]));
	test_it(($l[9],$l[10],$l[8]));
	test_it(($l[12],$l[13],$l[11]));
	test_it(($l[15],$l[16],$l[14]));
	test_it(($l[18],$l[19],$l[17]));
	test_it(($l[21],$l[22],$l[20]));
	test_it(($l[24],$l[25],$l[23]));
	chomp $d;
	chomp $p;
	chomp $c;
	chomp $o;
	chomp $f;
	chomp $g;
	chomp $s;
	print"$l[0],$d,$p,$c,$o,$f,$g,$s\n";
	$d ="unknown,100";
	$p ="unknown,100";
	$c ="unknown,100";
	$o ="unknown,100";
	$f ="unknown,100";
	$g ="unknown,100";
	$s ="unknown,100";

}

sub test_it {

	switch	($_[0]) {
		case "domain" {$d="$_[2],$_[1]";}
		case "phylum" {$p="$_[2],$_[1]";}
		case "class" {$c="$_[2],$_[1]";}
		case "order" {$o="$_[2],$_[1]";}
		case "family" {$f="$_[2],$_[1]";}
		case "genus" {$g="$_[2],$_[1]";}
		case "species" {$s="$_[2],$_[1]";}
	}
}
