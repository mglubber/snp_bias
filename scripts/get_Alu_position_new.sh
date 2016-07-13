#!/bin/bash
# get Alu positions from UCSC hg38 repeatMasker file. Include overlap of the ALU
# element with the database consensus sequence

# usage: "./get_Alu_positions_new.sh hg38.fa.out"
FILE="$1"
# Add headers to the top of the file for use in R. Respectively
# Chromosome(query) name, query start, query end, strand (+ or -), Repeat Family, Consensus start, end, missing
# Consensus positions start relative to consensus sequence (pos 1), end is final position in consensus, missing
# is the number of bp from end of match to the end of the consensus sequence (i.e. if the repeat match starts 
# from the 5 bp of a 215bp consensus sequence and ends at the 200th, Con_Start, Con_End, and Con_Miss will be 
# 5, 200, (15) respectively. If the repeate is in the reverse strand, the order of  Con_Start, Con_End, and 
# Con_Miss are reversed (i.e. (15), 200, 5)
echo "CHROM,START,END,STRAND,FAMILY,Con_Start,Con_End,Con_Miss" 
# grep the file for Alu repeats, replace all whitespace with tabs, change the negative strange from C to -,
# and remove unneccesary columns. 
grep 'Alu' "$FILE" | sed -r -e 's/^\s+//' -e 's/\s+/,/g' -e 's/,C,/,-,/'| cut -d "," -f 5-7,9-10,12-14
