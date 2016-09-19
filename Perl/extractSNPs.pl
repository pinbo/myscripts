#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( first );

# change here for different trait and different NAMs
my $trait = "TotalSN";
my $nam = "DD";

my @raw = glob("pvalues*$nam*");
my %pvalues;

open(OUT,">${trait}_${nam}.txt") or die "Cannot write to the file!\n";
print OUT "SNP";

foreach my $file (@raw){
	print OUT "\t$file";
	open (RAW, $file);
	chomp(my $header = <RAW>); # the first line
	my @ID = split(/\t/, $header);
	my $n = first {$ID[$_] eq $trait} 0 .. $#ID; # the column # of the trait
	while (<RAW>) {
		chomp(my $line = $_);
		my @col = split(/\t/, $line);
		my $key = $col[0];
		my $value = $col[$n];
		next if $value > 0.01;
		$pvalues{$key} .= "\t$value";
	}
}

print OUT "\n";


while ( my($key, $value) = each %pvalues ){
  my $count = () = $value =~ /\t/g; # count tab occurences
  if ($count > 1) {print OUT "$key$value\n";} # only SNPs significant in more than 1 environment were printed out
}

close OUT;




