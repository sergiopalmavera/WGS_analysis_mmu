#!/bin/bash

cat ../output/*_with_rec_rate.roh | grep -E "^# RG|^RG" > ../output/ROH_with_rec_rate_regions.table
