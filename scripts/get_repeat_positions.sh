#! /bin/bash
# get repeate positions from UCSC hg38 repeatMasker file. Include overlap of the repeat
# element with the database consensus sequence

# usage: "./get_repeat_positions.sh hg38.fa.out"
FILE="$1"
REPEAT_FILTER="$2"
# Add headers to the top of the file for use in R. Respectively
# Chromosome(query) name, query start, query end, strand (+ or -), Repeat Family, Consensus start, end, missing
# Consensus positions start relative to consensus sequence (pos 1), end is final position in consensus, missing
# is the number of bp from end of match to the end of the consensus sequence (i.e. if the repeat match starts 
# from the 5 bp of a 215bp consensus sequence and ends at the 200th, Con_Start, Con_End, and Con_Miss will be 
# 5, 200, (15) respectively. If the repeate is in the reverse strand, the order of  Con_Start, Con_End, and 
# Con_Miss are reversed (i.e. (15), 200, 5)
echo "CHROM,START,END,STRAND,FAMILY,Con_Start,Con_End,Con_Miss,ID" 
# grep the file for repeats matching the filter, replace all whitespace with tabs, change the negative strand
# from C to -, remove bracket from missing bp, and remove unneccesary columns. 
grep "$REPEAT_FILTER" "$FILE" | sed -r -e 's/^\s+//' -e 's/\s+/,/g' -e 's/,C,/,-,/' -e 's/\(([0-9]+)\)/\1/g' | cut -d "," -f 5-7,9-10,12-15
