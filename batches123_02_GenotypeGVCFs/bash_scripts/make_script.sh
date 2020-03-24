#!/bin/bash

chr=$1

new_script=$(echo GenotypeGVCFs_template.sh | sed "s/template/chr$chr"/) 

sed "s/chr=/chr=$chr/" GenotypeGVCFs_template.sh > $new_script

chmod +x $new_script

nohup $new_script &> ${new_script/.sh/.out} &
