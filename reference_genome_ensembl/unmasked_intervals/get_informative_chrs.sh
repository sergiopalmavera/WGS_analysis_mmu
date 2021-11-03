#!/bin/bash

# Extract chromosomes used for windowed Fst analysis followed by detection of RDDs.
# For that analysis only autosomes (chr 1 to 19) and the X chr (they are only females) were used.

grep -E "^[1-9]|^X" Mus_musculus.GRCm38.dna.primary_assembly.genome > Mus_musculus.GRCm38.dna.primary_assembly_chr1_to_19_plus_chrX.genome
