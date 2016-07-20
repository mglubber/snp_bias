library("GenomicRanges")
get_flank_distances <- function(query_range, subject_range, flank_distance, offset, upstream=FALSE){
    subject_positive <- subject_range[strand(subject_range)=="+"]
    subject_negative <- subject_range[strand(subject_range)=="-"]
    # find overlaps between snps and downstream region of positive strand alu
    # elements and snps, ignoring the 'strand' of the SNP.
    # flank (start=FALSE) defines region from end+1 to end+flankdistance.
    # +offset extends both sides of resulting range by offset. Combined gives
    # a interval starting offset bp before the 3' end and extending
    # flank_distance-offset bp downstream
    positive_overlaps <- findOverlaps(
        query_range, 
        (flank(subject_positive,width=flank_distance-(2*offset),start=upstream,))+offset,
        ignore.strand=TRUE)
    # repeat for negative strand alus. Genomic ranges automatically adjust flank
    # for different strands.
    negative_overlaps <- findOverlaps(
        query_range,
        (flank(subject_negative,width=flank_distance-(2*offset),start=upstream,))+offset,
        ignore.strand=TRUE)
    positive_distance <- cbind.data.frame(
        alu_family=subject_positive$alufamily[subjectHits(positive_overlaps)], 
        distance=(start(query_range)[queryHits(positive_overlaps)] - 
                  end(subject_positive)[subjectHits(positive_overlaps)]))
    negative_distance <- cbind.data.frame(
        alu_family=subject_negative$alufamily[subjectHits(negative_overlaps)], 
        distance=(start(subject_negative)[subjectHits(negative_overlaps)] - 
            end(query_range)[queryHits(negative_overlaps)]))
    #positive_distance <- positive_distance[positive_distance[,2]<(flank_distance+1) & positive_distance[,2]>-6,]
    #negative_distance <- negative_distance[negative_distance[,2]<(flank_distance+1) & positive_distance[,2]>-6,]
    distance_vs_family <- rbind.data.frame(positive_distance, negative_distance)
    result <- list(d_v=distance_vs_family,p_d=positive_distance,n_d=negative_distance)
    return(result)
}

result_list <- get_flank_distances(snp_ranges[snp_ranges$moltype=="genomic"&snp_ranges$class=="single",],alu_ranges, 200, 5, FALSE)
distance_vs_family <- result_list$d_v
positive_distance <- result_list$p_d
negative_distance <- result_list$n_d
