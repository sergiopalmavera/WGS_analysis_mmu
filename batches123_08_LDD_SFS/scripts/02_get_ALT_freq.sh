#!/bin/bash

for fl in ../output/*.frq
do
	echo $fl
	awk '{print $6}' $fl | awk -F ":" '{print $2}' | sed -n '2,$p' > ${fl/.frq/.ALTfrq}
	printf "\n"
done
