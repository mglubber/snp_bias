library("GenomicRanges")
source('scripts/overlap_distance_functions.R')
alu_file <- "data_int/ALU_locations.txt"
snp_file <- "data_int/allsnpsmodified.txt"

alu_table <- read.csv(alu_file,header=T)
cat("Created table of alu elements from file", alu_file, "\n")
cat(nrow(alu_table), "repeat elements found.\n")
# Filter for only repeats with unique IDs, to prevent biasing from repeats 
# duplicated through other mechanisms
alu_table <- alu_table[alu_table$ID %in% names(table(alu_table$ID))[table(alu_table$ID) <= 1],]
cat(nrow(alu_table), "unique repeat elements remaining.\n")
# flip Con_Start and Con_End for negative strand repeats, so that column names
# are correct
alu_table[alu_table$STRAND=="-",c(6,8)] <- alu_table[alu_table$STRAND=="-",c(8,6)]
cat("Flipped consensus sequence start positions and truncated base pairs for negative strand\n")
# Convert the table to a set of genomic ranges, with additional information
alu_ranges <- GRanges(
    seqname=alu_table$CHROM,
    ranges=IRanges(start=alu_table$START,end=alu_table$END),
    strand=alu_table$STRAND,
    alufamily=alu_table$FAMILY,
    conStart=alu_table$Con_Start,
    conEnd=alu_table$Con_End,
    conMiss=alu_table$Con_Miss)
cat("Created", length(alu_ranges), "repeat element genomic ranges\n")
# Read snp/other mutations locations from file
snp_table <- read.table(snp_file, sep='\t', header=T)
cat("Created table of", nrow(snp_table), "SNP/Mutation positions and types from", snp_file, ".\n")
# Create genomic ranges from snp/mutation tables
snp_ranges <- GRanges(
    seqname=snp_table$chrom,
    ranges=IRanges(start=snp_table$chromStart,end=snp_table$chromEnd),
    strand=snp_table$strand,
    name=snp_table$name,
    observed=snp_table$observed,
    moltype=snp_table$molType,
    class=snp_table$class)
cat("Created", length(snp_ranges), "SNP/Mutation genomic ranges.\n")

# Find distances between pairs of nearby alus and mutations using 
# get_flank_distances from the overlap_distance_functions.R script.
alu_all_pairs <- get_flank_distances(
    query_range=snp_ranges,subject_range=alu_ranges,
    flank_distance=200,offset=5,upstream=FALSE)

# Find distances for just SNP, a.k.a genomic singles in the snp_ranges
alu_snp_pairs <- get_flank_distance(
    query=snp_ranges[snp_ranges$class=="single"&snp_ranges$moltype=="genomic",
    subject_range=alu_ranges,flank_distance=200,offset=5,upstream=FALSE)

# Do the same for upstream SNPs
alu_snp_upstream <- get_flank_distance(
    query=snp_ranges[snp_ranges$class=="single"&snp_ranges$moltype=="genomic",
    subject_range=alu_ranges,flank_distance=200,offset=5,upstream=FALSE)
