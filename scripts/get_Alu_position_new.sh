#!/bin/bash
# get Alu positions from UCSC hg38 repeatMasker file. Include overlap of the ALU
# element with the database consensus sequence
FILE="$1"
echo "CHROM	START	END	SOMETHING	STRAND	FAMILY	E_MISS	E_START	E_END"
grep 'Alu' "$FILE" | sed -r -e 's/\s+/\t/g' -e 's/\tC\t/\t-\t/'| cut -f 6-11,13-15
