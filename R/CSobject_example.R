#' Create example CapSet object
#'
#' @return An object of class CapSet
#' @export
#'
#' @examples
#' cs <- exampleCSobject()
#'
exampleCSobject <- function() {
    # list of barcode IDs
     #idxlist <- c("CAAGTG", "AGATGC", "TGTGAG",
     #            "GGTTAC", "TTAGCC", "AGTCGA")
     # corresponding sample names
     #fnames <-  c("WTa", "WTb", "WTc",
     #            "MLEa", "MLEb", "MLEc")
    idxlist <- c("CAAGTG", "TTAGCC", "GTGGAA", "TGTGAG")
    fnames <- c("embryo1", "embryo2", "embryo3", "embryo4")
    # numbers
    demult_reads <- c(628, 778, 648, 291, 297, 311)
    num_mapped <- c(1232, 1512, 1266, 564, 574, 606)
    num_filtered <- c(569, 715, 584, 240, 266, 277)
    num_intss <- num_filtered

    dir <- system.file("extdata", package = "icetea")

    si <- data.frame(row.names = idxlist,
                samples = fnames,
                demult_reads = demult_reads,
                demult_R1 = file.path(dir, 'demulti_fastq', paste0(fnames, "_R1.fastq.gz")),
                demult_R2 = file.path(dir, 'demulti_fastq', paste0(fnames, "_R2.fastq.gz")),
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
     tss <- load(file.path(dir, 'testTSS.Rdata'))
     cs@tss_detected <- tss
     return(cs)
     }
