#!/bin/bash

zcat ../output/cohort.vcf.gz | grep -v "^#" | wc -l
zcat ../output/cohort_biallelicSNPs.vcf.gz | grep -v "^#" | wc -l


