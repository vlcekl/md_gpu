#!/usr/bin/perl -w
#
# File name:   makeco2.pl
# Date:        2009/06/03 01:56
# Author:      Lukas Vlcek
#
# Description: 
#

$ll = 1.163;
$qo = -0.3256;
$qc = 0.6512;

$n = 1000; 

$lx = 71;
$ly = 71;
$lz = 71;
$nx = 10;
$ny = 10;
$nz = 10;
$dx = $lx/$nx;
$dy = $ly/$ny;
$dz = $lz/$nz;
$xm = 0.0;
$ym = 0.0;
$zm = 0.0;

$xx = $ll/sqrt(3.0);

$imin = 0;

$i = $imin;
for $ix (1 .. $nx) {
    $rx = $ix*$dx + $xm;
    for $iy (1 .. $ny) {
        $ry = $iy*$dy + $ym;
        for $iz (1 .. $nz) {
            $rz = $iz*$dz + $zm;
            $i++;
            printf "O%24.14f%24.14f%24.14f\n", $rx+$xx, $ry+$xx, $rz+$xx;
            $i++;
            printf "C%24.14f%24.14f%24.14f\n", $rx, $ry, $rz;
            $i++;
            printf "O%24.14f%24.14f%24.14f\n", $rx-$xx, $ry-$xx, $rz-$xx;
        }
    }
}

