#!/bin/bash

cat ../output/*.roh | grep -E "^# RG|^RG" > ../output/ROH_regions.table
