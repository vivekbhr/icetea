#' Create example CapSet object
#'
#' @return An object of class CapSet
#' @export
#'
#' @examples
#' cs <- exampleCSobject(protocol = "MAPCap")
#'
exampleCSobject <- function(protocol = "MAPCap") {
    # list of barcode IDs
    idxlist <- c("CAAGTG", "TTAGCC", "GTGGAA", "TGTGAG")
    fnames <- c("embryo1", "embryo2", "embryo3", "embryo4")
    # numbers
    demult_reads <- c(60, 80, 180, 140)
    num_mapped <- c(56, 79, 178, 138)
    num_filtered <- c(47, 56, 145, 91)
    num_intss <- c(46, 53, 138, 88)

    dir <- system.file("extdata", package = "icetea")
    si <- data.frame(row.names = idxlist,
                samples = fnames,
                demult_reads = demult_reads,
                demult_R1 = file.path(dir, 'demult_fastq', paste0(fnames, "_R1.fastq.gz")),
                demult_R2 = file.path(dir, 'demult_fastq', paste0(fnames, "_R2.fastq.gz")),
                mapped_file = file.path(dir, 'bam', paste0(fnames, ".bam")),
                num_mapped = num_mapped,
                filtered_file = file.path(dir, 'filtered_bam', paste0(fnames, ".filtered.bam")),
                num_filtered = num_filtered,
                num_intss = num_intss
                )

     cs <- newCapSet(expMethod = 'MAPCap', fastqType = 'paired',
                    fastq_R1 = file.path(dir, 'mapcap_test_R1.fastq.gz'),
                    fastq_R2 = file.path(dir, 'mapcap_test_R2.fastq.gz'),
                    sampleInfo = si)
     # tss detected
     tss <- load(file.path(dir, "testTSS_grl.Rdata"))
     cs@tss_detected <- tss
     return(cs)
     }
