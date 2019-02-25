#!/usr/bin/perl -w
#
# File name:   mktop.pl
# Date:        2009/05/26 16:14
# Author:      Lukas Vlcek
#
# Description: 
#

$n2 = 4096;

for $i (1 .. $n2) {
    $j = ($i-1)*2 + 1;
    $k = ($i-1)*2 + 2;
    print "$j $k\n";
}

#for $i (1 .. $nwat) {
#    $j = ($i-1)*3 + 1;
#    $k = ($i-1)*3 + 2;
#    $l = ($i-1)*3 + 3;
#    print "$j $k $l\n";
#}




# end of mktop.pl 
