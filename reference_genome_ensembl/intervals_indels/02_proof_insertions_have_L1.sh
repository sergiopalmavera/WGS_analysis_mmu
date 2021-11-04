#!/bin/bash

# Are insertions always replacement of one base in the reference (it could be that two bases are replaced by insertion 50 bases!)
# Long story short: it's always one base being replace by an insertion

awk '{print $2 "\t" $3 "\t" $3-$2}' mus_musculus_ins.bed | awk '{print $3}' | uniq
