use strict;
use warnings;

my $term = join '|', map "\Q$_\E", split ' ', pop; # get the pattern, if more than one patterns, will be jointed with "|" for OR.
#print $term;
my $found = 0;

while (<>) {
    if (/^>/) {
        #$found = /$term/i ? 1 : 0;
        #print if $found;
        #next;
        $found = 0;
        if (/$term:(\d+):(\d+):/) {
        	#print "$1\n";
    		$found = ($1>159000000 and $1<165000000) ? 1 : 0;
    	}
    }

    print if $found;
}
