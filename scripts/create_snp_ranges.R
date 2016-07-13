library("GenomicRanges")
alu_table <- read.table("data/ALU_locations.txt",sep='\t',header=T)
alu_ranges <- GRanges(seqname=alu_table$CHROM,ranges=IRanges(start=alu_table$START,end=alu_table$END),strand=alu_table$STRAND,alufamily=alu_table$FAMILY)
alu_positive <- alu_ranges[strand(alu_ranges)=="+"]
alu_negative <- alu_ranges[strand(alu_ranges)=="-"]
snp_table <- read.table("genomicsproject/allsnpsmodified", sep='\t', header=T)
snp_ranges <- GRanges(seqname=snp_table$chrom,ranges=IRanges(start=snp_table$chromStart,end=snp_table$chromEnd),strand=snp_table$strand,name=snp_table$name)
