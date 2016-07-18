library("GenomicRanges")
alu_file <- "data_int/ALU_locations.txt"
snp_file <- "data_int/allsnpsmodified.txt"

alu_table <- read.csv(alu_file,header=T)
cat("Created table of alu elements from file", alu_file, "\n")
cat(nrow(alu_table), "repeat elements found.\n")
# flip Con_Start and Con_End for negative strand repeats, so that column names are correct
alu_table[alu_table$STRAND=="-",c(6,8)] <- alu_table[alu_table$STRAND=="-",c(8,6)]
cat("Flipped consensus sequence start positions and truncated base pairs for negative strand\n")
# Convert the table to a set of genomic ranges, with additional information
alu_ranges <- GRanges(seqname=alu_table$CHROM,ranges=IRanges(start=alu_table$START,end=alu_table$END),strand=alu_table$STRAND,alufamily=alu_table$FAMILY,conStart=alu_table$Con_Start,conEnd=alu_table$Con_End,conMiss=alu_table$Con_Miss)
cat("Created", length(alu_ranges), "repeat element genomic ranges\n")
# Split alu_ranges in positive and negative ranges, to simplify analysis of downstream/upstream SNPs
alu_positive <- alu_ranges[strand(alu_ranges)=="+"]
alu_negative <- alu_ranges[strand(alu_ranges)=="-"]
cat("Split positive and negative strand genomic ranges to simplify finding overlaps.\n Loading SNP/mutation data, this make take some time.\n")
# Read snp/other mutations locations from file
snp_table <- read.table(snp_file, sep='\t', header=T)
cat("Created table of", nrow(snp_table), "SNP/Mutation positions and types from", snp_file, ".\n")
# Create genomic ranges from snp/mutation tables
snp_ranges <- GRanges(seqname=snp_table$chrom,ranges=IRanges(start=snp_table$chromStart,end=snp_table$chromEnd),strand=snp_table$strand,name=snp_table$name,observed=snp_table$observed,moltype=snp_table$molType,class=snp_table$class)
cat("Created", length(snp_ranges), "SNP/Mutation genomic ranges.\n")
