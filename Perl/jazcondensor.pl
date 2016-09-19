#!/usr/bin/perl
#  User input | read file | write file 
#  Copyright 2012 Junli Zhang <junli@junli-ThinkPad-T510>
# To join all the Jaz data files into one file with the selected wavelength for each channel.
#  
use strict;
use warnings;

print "Please select on channel (0, 1, 2 ...):\n";
chomp(my $channel = <STDIN>);

my @raw = glob("Sp$channel*");
open(IN, "/home/junli/bin/SP$channel") or die "Cannot open the file SP0";
open(OUT,">SP${channel}_OUTPUT.txt") or die "Cannot write to the file";
my %seed;

foreach my $index (<IN>){
	chomp($index);
	$seed{$index} = $index;
}
my $num = keys(%seed);

print OUT "Index";

foreach my $file (@raw){
	print OUT "\t$file";
	my $count = 0;
	open (RAW, $file);
	while (<RAW>) {
		chomp(my $line = $_);
		my ($key, $value) = split(/\s+/, $line);
		if (exists $seed{$key}) {$seed{$key} .= "\t$value"; $count += 1;}
		last if $count == $num;		
	}
}
print OUT "\n";

foreach my $key (sort {$a <=> $b} keys(%seed)){print OUT "$seed{$key}\n"}

