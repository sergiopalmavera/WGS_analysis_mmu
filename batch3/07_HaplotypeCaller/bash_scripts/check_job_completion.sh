#!/bin/bash

# This script checks the completion of haplotypecaller by searching the final tool's message "org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done" while manking sure that the the time of completion is reasonable (thousands of minutes!)


for f in ./*.out
do
	echo "Batch $f"
	grep "org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done." $f
	printf "\n"
done > check_job_completion.txt
