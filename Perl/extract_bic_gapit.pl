#!/usr/bin/env perl
#  temp.pl
#  
#  Copyright 2014 Junli Zhang <zhjl86@gmail.com>
use warnings;
use strict;

my @files = glob("*.BIC.Model.Selection.Results.csv");
open(OUT,">Conbined_BIC.csv") or die "Cannot write to the file";
print OUT "Trait,Number of CV,BIC,Loglik\n";
foreach my $file (@files){
	my $trait = substr $file, 7, -32;
	open(IN, $file);
	my $header = <IN>; #skip header
	#only print the first line and last line	
	my $first = <IN>;
	my $last;
	while (<IN>) { $last = $_ }
	print OUT "$trait,$first";
	print OUT "$trait,$last";
}
close OUT;
