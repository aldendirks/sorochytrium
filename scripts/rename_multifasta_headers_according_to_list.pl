# Source: https://www.biostars.org/p/103089/
# Usage: perl script.pl textFile fastaFile [>outFile]

use strict;
use warnings;

my @arr;

while (<>) {
    chomp;
    push @arr, $_ if length;
    last if eof;
}

while (<>) {
    print /^>/ ? shift(@arr) . "\n" : $_;
}