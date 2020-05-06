#!/bin/bash

chr=$1

template=concatenate_sort_validate_template.sh
script=${template/template/chr$chr}

sed 's/chr=/chr='$1'/' $template > $script

chmod +x $script 

nohup  $script &> ${script/.sh/.out} &

sleep 20

ps -xjf | head -1 > tmp1
ps -xjf | grep $script | grep bash > tmp2

cat tmp1 tmp2 > ${script/.sh/.pid}

rm tmp1 tmp2

