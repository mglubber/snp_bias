#!/bin/bash
# get Alu positions from UCSC hg38 repeatMasker file. Include overlap of the ALU
# element with the database consensus sequence
FILE="$1"
echo "CHROM	START	END	SOMETHING	STRAND	FAMILY	E_MISS	E_START	E_END" > alu_positions.tab
tail -n+4 "$FILE" | sed -n '/SINE\/Alu/p' | sed -e 's/^\s\+//' -e 's/\s\+/\t/g' \
	-e 's/\sC\s/\t-\t/' | cut -f 5-10,12-14 >> alu_positions.tab

