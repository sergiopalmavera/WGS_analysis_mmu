#!/bin/bash

chr=$1
from=$2
to=$3
batch_size=$4

nm=$(echo consolidate_by_intervals_chr00_fromXtoY_template.sh | sed 's/chr00/chr'$chr'/ ; s/fromX/from'$from'/ ; s/toY/to'$to'/ ; s/_template//')
echo "# New script name: $nm"

sed 's/chr=/chr='$chr'/ ; s/from=/from='$from'/ ; s/to=/to='$to'/ ; s/batch_size=/batch_size='$batch_size'/' ./consolidate_by_intervals_chr00_fromXtoY_template.sh > $nm

chmod +x $nm

./$nm &> ./${nm/.sh/.out}
