#! /bin/bash
# get summary of 'snp' types (e.g. genomic_single, genomic_insertion, etc) from
# file and sort from most to least common
FILE = $1
awk '{print $11"_"$12}' "$FILE" | sort | uniq -c | sort -nr
