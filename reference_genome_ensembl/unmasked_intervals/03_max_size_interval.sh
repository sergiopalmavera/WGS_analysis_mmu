#!/bin/bash

awk '{print $3 - $2}' intervals_unmasked_excl_indels.bed | sort -hr | uniq | head
