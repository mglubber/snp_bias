library("GenomicRanges")
flank_distance <- 200
positive_genomic_overlaps <- findOverlaps(snp_ranges,(flank(alu_positive,width=flank_distance,start=FALSE)),type="within",ignore.strand=TRUE)
negative_genomic_overlaps <- findOverlaps(snp_ranges,(flank(alu_negative,width=flank_distance,start=FALSE)),type="within",ignore.strand=TRUE)
positive_distance <- cbind.data.frame(alu_family=alu_positive$alufamily[subjectHits(positive_genomic_overlaps)], distance=(start(snp_ranges)[queryHits(positive_genomic_overlaps)] - end(alu_positive)[subjectHits(positive_genomic_overlaps)]))
negative_distance <- cbind.data.frame(alu_family=alu_negative$alufamily[subjectHits(negative_genomic_overlaps)], distance=(start(alu_negative)[subjectHits(negative_genomic_overlaps)] - end(snp_ranges)[queryHits(negative_genomic_overlaps)]))
positive_distance <- positive_distance[positive_distance[,2]<(flank_distance+1),]
negative_distance <- negative_distance[negative_distance[,2]<(flank_distance+1),]
distance_vs_family <- rbind.data.frame(positive_distance, negative_distance)

