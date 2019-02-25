#!/usr/bin/perl -w
#
# File name:   makeco2.pl
# Date:        2009/06/03 01:56
# Author:      Lukas Vlcek
#
# Description: 
#

$ll = 1.0976;
$qo = -0.3256;
$qc = 0.6512;

$n = 4096; 
$nat = 3*$n;

$lx = 170;
$ly = 170;
$lz = 170;
$nx = 16;
$ny = 16;
$nz = 16;
$dx = $lx/$nx;
$dy = $ly/$ny;
$dz = $lz/$nz;
$xm = -$lx/2;
$ym = -$ly/2;
$zm = -$lz/2;

$xx = $ll/sqrt(3.0);

$imin = 0;

print "$nat 1\n";
print "$lx $ly $lz\n";

$i = $imin;
for $ix (1 .. $nx) {
    $rx = $ix*$dx + $xm;
    for $iy (1 .. $ny) {
        $ry = $iy*$dy + $ym;
        for $iz (1 .. $nz) {
            $rz = $iz*$dz + $zm;

            $i++;
            print "$i 2 $qo\n";
            $ax = $rx - $xx;
            $ay = $ry - $xx;
            $az = $rz - $xx;
            print "$ax $ay $az\n";
            print "0.0 0.0 0.0\n";

            $i++;
            print "$i 1 $qc\n";
            print "$rx $ry $rz\n";
            print "0.0 0.0 0.0\n";

            $i++;
            print "$i 2 $qo\n";
            $ax = $rx + $xx;
            $ay = $ry + $xx;
            $az = $rz + $xx;
            print "$ax $ay $az\n";
            print "0.0 0.0 0.0\n";
        }
    }
}

# BONDS


print "# Bonds\n";

$ii = 0;
for $i (1 .. $n) {
    $j = ($i-1)*3 + 1;
    $k = ($i-1)*3 + 2;
    $l = ($i-1)*3 + 3;
    $ii++;
    print "$ii 1 $j $k\n";
    $ii++;
    print "$ii 1 $k $l\n";
}

print "# Angles\n";

$ii = 0;
for $i (1 .. $n) {
    $j = ($i-1)*3 + 1;
    $k = ($i-1)*3 + 2;
    $l = ($i-1)*3 + 3;
    $ii++;
    print "$ii 1 $j $k $l\n";
}
