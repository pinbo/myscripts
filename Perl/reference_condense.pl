#!/usr/bin/perl
#  Usage: name.pl reference_prefix (Sdark0 OR Sdark1 etc)
#  Copyright 2014 Junli Zhang <zhjl86@gmail.com>
#  To join raw reference files.
#  
use strict;
use warnings;

#print "Please select on channel (0, 1, 2 ...):\n";
#chomp(my $channel = <STDIN>);

# Make sure the correct number of command-line args
#my $num_args = $#ARGV + 1;
if (@ARGV != 1) { # only 1 argument
    print "\nUsage: name.pl reference_prefix (Sdark0 OR Sdark1 etc)\n\n";
    exit;
}

my $channel=$ARGV[0]; # Sdark0 / Swhite0 etc

my @raw = glob("$channel*");
open(OUT,">Condensed_${channel}.txt") or die "Cannot write to the file!\n";
my %csr;

print OUT "Index\t", join("\t", @raw), "\n";

foreach my $file (@raw){
	open (RAW, $file);
	while (<RAW>) {
		chomp(my $line = $_);
		my ($key, $value) = split(/\s+/, $line);
		$csr{$key} .= "\t$value";
	}
}

foreach my $key (sort {$a <=> $b} keys(%csr)){print OUT "$key$csr{$key}\n"};

close OUT;

