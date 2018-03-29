#' Create example CapSet object
#'
#' @param expMethod Which experiment method to use (options : "RAMPAGE", "MAPCap")
#' @return An object of class CapSet
#' @export
#'
#' @examples
#' cs <- exampleCSobject(expMethod = "MAPCap")
#'
exampleCSobject <- function(expMethod = "MAPCap") {
    # list of barcode IDs
    idxlist <- c("CAAGTG", "TTAGCC", "GTGGAA", "TGTGAG")
    fnames <- c("embryo1", "embryo2", "embryo3", "embryo4")
    # numbers
    demult_reads <- c(60, 80, 180, 140)
    num_mapped <- c(56, 79, 178, 138)
    num_filtered <- c(47, 56, 145, 91)
    # reads in tss are the expected reads when TSScalling is done with 2-fold enrichment
    num_intss <- c(23, 32, 41, 40)

    dir <- system.file("extdata", package = "icetea")
    cs <- newCapSet(
        expMethod = 'MAPCap',
        fastq_R1 = file.path(dir, 'mapcap_test_R1.fastq.gz'),
        fastq_R2 = file.path(dir, 'mapcap_test_R2.fastq.gz'),
        idxList = idxlist,
        sampleNames = fnames,
        demult_R1 = file.path(dir, 'demult_fastq', paste0(fnames, "_R1.fastq.gz")),
        demult_R2 = file.path(dir, 'demult_fastq', paste0(fnames, "_R2.fastq.gz")),
        mapped_file = file.path(dir, 'bam', paste0(fnames, ".bam")),
        filtered_file = file.path(dir, 'filtered_bam', paste0(fnames, ".filtered.bam"))
    )
    si <- sampleInfo(cs)
    si$num_intss <- num_intss
    sampleInfo(cs) <- si
    # tss detected
    grl <- readRDS(file.path(dir, "testTSS_grl.rds"))
    cs@tss_detected <- grl
    return(cs)
}
