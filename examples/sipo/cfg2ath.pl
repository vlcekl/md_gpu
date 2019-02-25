#!/usr/bin/perl -w
#
# File name:   cfg2ath.pl
# Date:        2009/06/28 18:50
# Author:      Lukas Vlcek
#
# Description: 
#

$nw = 311;
$nli = 5;
$ncl = 5;
$nsi = 1072;
$no = 2117;
$noh = 54;

$n = 3*$nw + $nli + $ncl + $nsi + $no + 2*$noh;

$l=<> for (1 .. 2);
$l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
$lx = $1; $ly = $2; $lz = $3;
$l=<> for (1 .. 6);

print "$n 1\n";
print "$lx $ly $lz\n";

$i = 0;
for (1 .. $nw) {
    $l=<> for (1 .. 2);
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $x = $1; $y = $2; $z = $3;
    $l=<>;
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $vx = $1; $vy = $2; $vz = $3;
    $l=<> for (1 .. 3);
    $i++;
    print "$i 2 0.41\n";
    print "$x $y $z\n";
    print "$vx $vy $vz\n";

    $l=<> for (1 .. 2);
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $x = $1; $y = $2; $z = $3;
    $l=<>;
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $vx = $1; $vy = $2; $vz = $3;
    $l=<> for (1 .. 3);
    $i++;
    print "$i 1 -0.82\n";
    print "$x $y $z\n";
    print "$vx $vy $vz\n";

    $l=<> for (1 .. 2);
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $x = $1; $y = $2; $z = $3;
    $l=<>;
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $vx = $1; $vy = $2; $vz = $3;
    $l=<> for (1 .. 3);
    $i++;
    print "$i 2 0.41\n";
    print "$x $y $z\n";
    print "$vx $vy $vz\n";
}

$l=<> for (1 .. 2);
    
for (1 .. $nli) {
    $l=<> for (1 .. 2);
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $x = $1; $y = $2; $z = $3;
    $l=<>;
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $vx = $1; $vy = $2; $vz = $3;
    $l=<> for (1 .. 3);
    $i++;
    print "$i 3 1.0\n";
    print "$x $y $z\n";
    print "$vx $vy $vz\n";
}
$l=<> for (1 .. 2);

for (1 .. $ncl) {
    $l=<> for (1 .. 2);
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $x = $1; $y = $2; $z = $3;
    $l=<>;
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $vx = $1; $vy = $2; $vz = $3;
    $l=<> for (1 .. 3);
    $i++;
    print "$i 4 -1.0\n";
    print "$x $y $z\n";
    print "$vx $vy $vz\n";
}
$l=<> for (1 .. 2);

for (1 .. $nsi) {
    $l=<> for (1 .. 2);
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $x = $1; $y = $2; $z = $3;
    $l=<>;
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $vx = $1; $vy = $2; $vz = $3;
    $l=<> for (1 .. 3);
    $i++;
    print "$i 5 2.1\n";
    print "$x $y $z\n";
    print "$vx $vy $vz\n";
}
$l=<> for (1 .. 2);

for (1 .. $no) {
    $l=<> for (1 .. 2);
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $x = $1; $y = $2; $z = $3;
    $l=<>;
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $vx = $1; $vy = $2; $vz = $3;
    $l=<> for (1 .. 3);
    $i++;
    print "$i 1 -1.05\n";
    print "$x $y $z\n";
    print "$vx $vy $vz\n";
}
$l=<> for (1 .. 2);

for (1 .. $noh) {
    $l=<> for (1 .. 2);
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $x = $1; $y = $2; $z = $3;
    $l=<>;
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $vx = $1; $vy = $2; $vz = $3;
    $l=<> for (1 .. 3);
    $i++;
    print "$i 1 -0.95\n";
    print "$x $y $z\n";
    print "$vx $vy $vz\n";

    $l=<> for (1 .. 2);
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $x = $1; $y = $2; $z = $3;
    $l=<>;
    $l=~/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/;
    $vx = $1; $vy = $2; $vz = $3;
    $l=<> for (1 .. 3);
    $i++;
    print "$i 2 0.425\n";
    print "$x $y $z\n";
    print "$vx $vy $vz\n";
}
$l=<> for (1 .. 4);

$nb = $nw*2 + $noh*2;
print "\nBonds: $nb\n\n";

$i = 0;
for $ii (1 .. $nw) {
    $j = ($ii-1)*3 + 1;
    $k = ($ii-1)*3 + 2;
    $l = ($ii-1)*3 + 3;
    $i++;
    print "$i 1 $j $k\n";
    $i++;
    print "$i 1 $k $l\n";
}

$istart = 4132;
$nstart = 1961;
for $ii (1 .. $noh) {
    $j = $nstart + $ii;
    $k = $istart + ($ii-1)*2 + 1;
    $l = $istart + ($ii-1)*2 + 2;
    $i++;
    print "$i 2 $j $k\n";
    $i++;
    print "$i 1 $k $l\n";
}

$nb = $nw + $noh;
print "\nAngles: $nb\n\n";
$i = 0;
for $ii (1 .. $nw) {
    $j = ($ii-1)*3 + 1;
    $k = ($ii-1)*3 + 2;
    $l = ($ii-1)*3 + 3;
    $i++;
    print "$i 1 $j $k $l\n";
}

$istart = 4132;
$nstart = 1961;
for $ii (1 .. $noh) {
    $j = $nstart + $ii;
    $k = $istart + ($ii-1)*2 + 1;
    $l = $istart + ($ii-1)*2 + 2;
    $i++;
    print "$i 2 $j $k $l\n";
}

# end of cfg2ath.pl 
