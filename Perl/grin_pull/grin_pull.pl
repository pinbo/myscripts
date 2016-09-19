#!/usr/bin/perl 
#http://www.ars-grin.gov/cgi-bin/npgs/acc/search.pl?accid=CItr+15471
use strict;
use warnings;

use LWP::Simple; 
sub grin{
	my $acc = shift @_;
	my $web = 'http://www.ars-grin.gov/cgi-bin/npgs/acc/search.pl?accid=';
	my $page = get($web . $acc);
	die "Couldn't get page" unless $page;
	my $receive_date = "";
	my $status = "";
	my $release = "";
	#$page =~ m{NPGS received: (\d?\d?-?[A-Z][a-z][a-z]-\d\d\d\d)\.}
	if ($page =~ m{NPGS received: ([^\.]*)\.}) {#anything but . OR until the first "."
		$receive_date = $1;
	} #NPGS received: 17-Jan-1962.
	if ($page =~ m{Improvement status: ([^\.]*)\.}) {
		$status = $1;
	} #Improvement status: Breeding material.
	if ($page =~ m{Released: ([^\.]*)\.}) {
		$release = $1;
	}#Released: 1973
	print "Completed Accessoin $acc\n";
	return "$receive_date\t$status\t$release";
}

# read input file
my $infile = 'cultivars.txt';
my $outfile = "release_date_cultivars.txt";
open(IN, $infile) or die "Could not open $infile: $!";
open(OUT,">$outfile") or die "Cannot write to the file";

while(<IN>)  {   
    chomp(my $index = $_);
	#print "$index\n";
	my $date = &grin($index);
	print OUT "$index\t$date\n";
}

close IN;
close OUT;