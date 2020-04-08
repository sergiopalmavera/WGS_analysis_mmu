#!/bin/bash

BEDTOOLS=/home/fb4/palma-vera/FBN_HOME/Tools/bedtools_version_2.29.2 #absolute path, modify accordingly

max_size=5000000

$BEDTOOLS/bedtools.static.binary makewindows -b ./intervals_unmasked.bed -w $max_size > intervals_unmasked_windowed_5Mb.bed
