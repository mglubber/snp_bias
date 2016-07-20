#! /bin/sh
# Find number of repeats with the same ID by family
 
$FILE=$1
$FILTER=$2

echo "#_of_alus #_of_fragments #family"
grep "$FILTER" "$FILE" | sed -r -e 's/^\s+//' -e 's/\s+/ /g' |\
awk '{ print $15 " " $10 }' | sort -k 1,1 -n | uniq -c |\
awk '{ print $1 " " $3 }' | sort -nr | uniq -c
