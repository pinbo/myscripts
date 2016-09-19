#!/usr/bin/perl 
#https://npgsweb.ars-grin.gov/gringlobal/accessiondetail.aspx?id=1154874

# try to use HTML::TableExtract
# because peopal do not recommend using regex to parse website
# Fedora install: sudo dnf install perl-HTML-TableExtract

# need to do more work

use strict;
use warnings;
use LWP::UserAgent; # for https
use HTML::TableExtract;


#my $ua = LWP::UserAgent->new;
#my $resp = $ua->get( 'https://npgsweb.ars-grin.gov/gringlobal/accessiondetail.aspx?id=1699131' );
#my $page = $resp->decoded_content;

## below code first help me find the table coordinate I want
#my $te = HTML::TableExtract->new;
#$te->parse($page);

#foreach my $ts ($te->tables) {
   #print "Table (", join(',', $ts->coords), "):\n";
   #foreach my $row ($ts->rows) {
      #print join(',', @$row), "\n";
   #}
#}

## then I can use this code to extract the table I need
#my $te = HTML::TableExtract->new( depth => 0, count => 8 );
#my $te = HTML::TableExtract->new( attribs => { class=>'grid horiz' } );
#$te->parse($page);
#foreach my $row ($te->rows) {
   #print join(',', @$row), "\n";
#}

# function
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
	my $te = HTML::TableExtract->new( attribs => { class=>'grid horiz' } );
	$te->parse($page);
	my @infor = ();
	foreach my $row ($te->rows) {push @infor, ${@$row}[1]};
	
	print "Completed Accession $acc\n";
	return @infor;
}

# read input file
my $infile = 'rust875ID_test.txt';
my $outfile = "rust875_GRIN_test.txt";
open(IN, $infile) or die "Could not open $infile: $!";
open(OUT,">$outfile") or die "Cannot write to the file";

while(<IN>)  { 
    chomp(my $index = $_);
	my $data = join "\t", &grin($index);
	print OUT "$index\t$data\n";
}

close IN;
close OUT;

