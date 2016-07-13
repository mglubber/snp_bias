#! /bin/bash
# find the counts of specific mutations (A/T, -/C, etc) from the allsnp file,
# sorted from most common to least
FILE=$1
awk '{print $10" "$11"_"$12}' $1 | sort | uniq -c | sort -nr
