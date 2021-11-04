REF=Mus_musculus.GRCm38.dna.primary_assembly.fa
echo "# Total Number of Sequences in reference $REF:"
grep "^>" $REF | wc -l
echo "# Done"
