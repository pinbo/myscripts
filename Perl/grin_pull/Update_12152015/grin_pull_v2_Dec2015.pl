#!/usr/bin/perl 
#https://npgsweb.ars-grin.gov/gringlobal/accessiondetail.aspx?id=1154874

# I need to install perl-LWP-Protocol-https in fedora
# sudo dnf install perl-LWP-Protocol-https
# it also installed perl-Mozilla-CA automatically

# check more examples on
# http://perlmaven.com/lwp-useragent-and-basic-authentication

use strict;
use warnings;
use LWP::UserAgent; # for https


sub grin{
	my $acc = shift @_;
	my $web = 'https://npgsweb.ars-grin.gov/gringlobal/accessiondetail.aspx?id=';
	my $ua = LWP::UserAgent->new;
	my $resp = $ua->get($web . $acc);
	my $page ="";
	if ($resp->is_success) {
        $page = $resp->decoded_content;
    } else {
        print $resp->status_line . "\n";
    }
	
	#print "Download of webpage done!\n";
	my $pedigree = "";
	my $status = "";
	if ($page =~ m{<th>Pedigree:</th>\s+<td>(.*?)</td>}) {
		$pedigree = $1;
		print "Pedigree: $pedigree!\n";
	} 
	if ($page =~ m{<th>Improvement status:</th>\s+<td>(.*?)</td>}) {
		$status = $1;
		print "Stutus: $status!\n";
	}
	print "Completed Accession $acc\n";
	return "$pedigree\t$status";
}

# read input file
my $infile = 'rust875ID.txt';
my $outfile = "rust875_pedigree.txt";
open(IN, $infile) or die "Could not open $infile: $!";
open(OUT,">$outfile") or die "Cannot write to the file";

while(<IN>)  { 
    chomp(my $index = $_);
	my $data = &grin($index);
	print OUT "$index\t$data\n";
}

close IN;
close OUT;
