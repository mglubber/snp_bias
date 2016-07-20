library("GenomicRanges")
flank_distance <- 50
positive_genomic_overlaps <- findOverlaps(snp_ranges,(flank(alu_positive,width=flank_distance-50,start=TRUE))+50L,type="within",ignore.strand=TRUE)
negative_genomic_overlaps <- findOverlaps(snp_ranges,(flank(alu_negative,width=flank_distance-50,start=TRUE))+50L,type="within",ignore.strand=TRUE)
positive_distance <- cbind.data.frame(alu_family=alu_positive$alufamily[subjectHits(positive_genomic_overlaps)], distance=(start(alu_positive)[subjectHits(positive_genomic_overlaps)] - end(snp_ranges)[queryHits(positive_genomic_overlaps)]))
negative_distance <- cbind.data.frame(alu_family=alu_negative$alufamily[subjectHits(negative_genomic_overlaps)], distance=(start(snp_ranges)[queryHits(negative_genomic_overlaps)] - end(alu_negative)[subjectHits(negative_genomic_overlaps)]))
distance_vs_family <- rbind.data.frame(positive_distance, negative_distance)

