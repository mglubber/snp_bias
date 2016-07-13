library("GenomicRanges")
flank_distance <- 200
alu_table <- read.table("data/ALU_locations.txt",sep='\t',header=T)
alu_ranges <- GRanges(seqname=alu_table$CHROM,ranges=IRanges(start=alu_table$START,end=alu_table$END),strand=alu_table$STRAND,alufamily=alu_table$FAMILY)
alu_positive <- alu_ranges[strand(alu_ranges)=="+"]
alu_negative <- alu_ranges[strand(alu_ranges)=="-"]
snp_table <- read.table("genomicsproject/allsnpsmodified", sep='\t', header=T)
snp_ranges <- GRanges(seqname=snp_table$chrom,ranges=IRanges(start=snp_table$chromStart,end=snp_table$chromEnd),strand=snp_table$strand,name=snp_table$name)
positive_genomic_overlaps <- findOverlaps(snp_ranges,(flank(alu_positive,width=flank_distance,start=FALSE)),type="within",ignore.strand=TRUE)
negative_genomic_overlaps <- findOverlaps(snp_ranges,(flank(alu_negative,width=flank_distance,start=FALSE)),type="within",ignore.strand=TRUE)
positive_distance <- cbind.data.frame(alu_family=alu_positive$alufamily[subjectHits(positive_genomic_overlaps)], distance=(start(snp_ranges)[queryHits(positive_genomic_overlaps)] - end(alu_positive)[subjectHits(positive_genomic_overlaps)]))
negative_distance <- cbind.data.frame(alu_family=alu_negative$alufamily[subjectHits(negative_genomic_overlaps)], distance=(start(alu_negative)[subjectHits(negative_genomic_overlaps)] - end(snp_ranges)[queryHits(negative_genomic_overlaps)]))
positive_distance <- positive_distance[positive_distance[,2]<(flank_distance+1),]
negative_distance <- negative_distance[negative_distance[,2]<(flank_distance+1),]
distance_vs_family <- rbind.data.frame(positive_distance, negative_distance)

