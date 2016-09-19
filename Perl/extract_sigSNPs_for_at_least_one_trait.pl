#! /usr/bin/perl

# Select SNPs that are signficant for at least one trait;
# for NAMs only

use warnings;
use strict; 
use List::Util qw( min );

# Make sure the correct number of command-line args
#if (@ARGV != 1) { # only 1 argument
    #print "\nUsage: program.pl input_file \n\n";
    #exit;
#}
#my $infile = $ARGV[0]; # input file

# function to extract signficant SNPs for each pvalue file
sub extract{
	my $infile = shift;
	open(OUT,">SigSNPs_$infile") or die "Cannot write to the file!\n";
	open (RAW, $infile);
	my $header = <RAW>; 
	print OUT $header;

	while(<RAW>){
		chomp(my $line = $_);
		my @pp = split(/\t/, $line);
		my $snp = shift @pp;
		my $min = min @pp;
		next if $min > 0.01;
		print OUT "$line\n";
	}
	close OUT;
}


my @raw = glob("pvalues*"); #pvalues of each NAM
foreach my $file (@raw) {&extract($file);}
