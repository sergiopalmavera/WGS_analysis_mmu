# Take the gtf gene set file and extract rows for genes, then ouput only postion, strand and gene name and the gene biotype
# To make it easier to process in R
awk '{ if ($3 == "gene") { print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $10 "\t" $14 "\t" $18} }' Mus_musculus.GRCm38.93.gtf | sed 's/\;//g ; s/"//g ; s/"//g' > Mus_musculus.GRCm38.93_gene.gtf
