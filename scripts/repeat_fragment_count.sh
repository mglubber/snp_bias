#! /bin/bash
# Usage: repeat_fragment_count.sh repeatfile repeatfilter
# Count the number of repeats with a given ID number from a repeatMasker output
# file. Assumes that repeats on a different chromosomes will never have the 
# same ID number.
FILE="$1"
REPEAT_FILTER="$2"
# Create header
echo "#_of_$REPEAT_FILTER #_of_fragments"
# search repeat file for repeats matching given filter
grep "$REPEAT_FILTER" "$FILE" |\
# remove/modify spacing for easier reading with awk
sed -r -e 's/^\s+//' -e 's/\s+/\t/g' |\
# split columns based on tabs, grab repeat IDs # and sort them numerically
awk 'BEGIN { FS = "\t" } ; { print $15 }' | sort -n
# count and sort by occurences of each ID number, then summarize number of occurences
uniq -c | sort -nr | awk '{ print $1 }' | uniq -c
