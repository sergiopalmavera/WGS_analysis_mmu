Pay no atention to the message "## something went wrong while bqsr" in the outputs.

This was suppose to be a message to warn about the final bam not found.

But I wrote the expression wrong:

- Instead of: "if [ -e ../output/${bam/.bam/.RG.dedup.bqsr.bam} ]"

- I wrote: "if [ -e $RES/${bam/.bam/.RG.dedup.bqsr.bam} ]", in which the directory $RES has not been specified.

This is a bug, but not a serious one....pay no attention to it.

Jobs were run on:
- 1 to 60 #linu2
- 61 to 70 #linu3
- 71 to 80 #linu5
- 81 to 90 #linu6

