grep "Alu" hg38.fa.out | sed -r s/s+/t/g | cut --complement -f2,12 | column -t > ALU_data_extended.txt
