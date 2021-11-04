# Before doing anything, make sure files have been produced properly:

echo "# Are entries in unmasked minus indels unique?"
wc -l ../unmasked_intervals/intervals_unmasked_minus_indels.bed
uniq ../unmasked_intervals/intervals_unmasked_minus_indels.bed | wc -l
printf "\n"

echo "# Are the entries in indels unique?"
wc -l ../intervals_indels/intervals_indels_sorted_uniq.bed 
uniq ../intervals_indels/intervals_indels_sorted_uniq.bed  | wc -l


