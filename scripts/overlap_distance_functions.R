library("GenomicRanges")

calculate_distance <- function(hits_df,upstream=FALSE){
# function that calculates the distance between two overlapping ranges.
# Can't use distance() function from GenomicRanges, because that package 
# considers overlapping ranges to have distance of 0. In this case, we
# want negative distances if the query 5' occurs before the subject 3'
# (or vice versa for upstream distances).

# Note: this function is very slow compared to the previous method of
# splitting into separate + and - strand lists and calculating distances
# separately. 

    attach(hits_df)
    if (upstream==TRUE) {
        hits_df$distance <- cbind(
            distance=ifelse(s_strand=="+",s_start-q_end,q_start-s_end))
        }
    else {
        hits_df$distance <- cbind(
            distance=ifelse(s_strand=="+",q_start-s_end,s_start-q_end))
        }
    
    return(hits_df)
    }

get_flank_distances <- function(query_range, subject_range, flank_distance, offset, upstream=FALSE){

    # find overlaps between snps and downstream region of repeat elements.
    # ignoring the 'strand' of the SNP.
    # flank (start=FALSE) defines region from end+1 to end+flankdistance.
    # +offset extends both sides of resulting range by offset. Combined gives
    # a interval starting offset bp before the 3' end and extending
    # flank_distance bp downstream

    flanking_overlaps <- findOverlaps(
        query_range, 
        (flank(subject_range,width=flank_distance-(offset),start=upstream,))+offset,
        ignore.strand=TRUE)

    # Get index positions of the repeat elements and snp in the genomic range
    # lists that overlap in the flank regions.
    subject_range_hits <- subjectHits(flanking_overlaps)
    query_range_hits <- queryHits(flanking_overlaps)

    # Create a data frame with relevant data from each pair of ranges that overlap.
    hits_df <- cbind.data.frame(
        s_family=subject_range$alufamily[subject_range_hits],
        s_strand=strand(subject_range[subject_range_hits]),
        s_start=start(subject_range[subject_range_hits]),
        s_end=end(subject_range[subject_range_hits]),
        q_name=query_range$name[query_range_hits],
        q_start=start(query_range[query_range_hits]),
        q_end=end(query_range[query_range_hits]),
        q_moltype=query_range$moltype[query_range_hits],
        q_class=query_range$class[query_range_hits],
        q_observed=query_range$observed[query_range_hits])
    
    # Calculate the distance and add to data frame. 
    hits_df <- calculate_distance(hits_df,upstream)

    # return data frame
    return(hits_df)
}
