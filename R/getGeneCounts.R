#' Get gene-level counts from TSS data
#'
#' @rdname getGeneCounts
#' @param CSobject The \code{\link{CapSet}} object to use.
#' @param transcriptGRL A GRangesList object containing transcripts, created using transcriptsBy(txdb)
#' @param regionAroundTSS integer, indicating how many bases downstream of TSS to count
#' @param outfile character. Tab-separated output file name (if required)
#' @param ncores integer. No. of cores/threads to use
#'
#' @return data.frame with gene-level counts for all genes in the txdb object
#'
#' @importFrom GenomicFeatures transcriptsBy
#'
#' @export
#' @examples
#'
#'  # load a txdb object
#'  library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#'  seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
#'
#'  # get transcripts by gene (only X chromsome, for simplicity)
#'  seqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "X"
#'  dm6trans <- transcriptsBy(TxDb.Dmelanogaster.UCSC.dm6.ensGene, "gene")
#'
#'  # load a CapSet object
#'  cs <- exampleCSobject()
#'  # get gene counts, counting reads around 500 bp of the TSS
#'  gcounts <- getGeneCounts(cs, dm6trans)
#'

setMethod("getGeneCounts",
        signature = "CapSet",
        function(CSobject,
                transcriptGRL,
                regionAroundTSS,
                outfile,
                ncores) {

        # get XX bp region around transcripts to count the reads
        transcriptGR <-
            GenomicRanges::resize(unlist(transcriptGRL), regionAroundTSS, fix = "start")
        # get bamfiles
        si <- sampleInfo(CSobject)
        if (all(is.na(si$filtered_file))) {
            warning("Filtered files not found under sampleInfo(CSobject). Using mapped files")
            bam.files <- si$mapped_file
        } else {
            bam.files <- si$filtered_file
        }

        # count reads
        bpParams <- getMCparams(ncores)
        tsscounts <- GenomicAlignments::summarizeOverlaps(transcriptGR,
                                                          reads = Rsamtools::BamFileList(bam.files),
                                                          singleEnd = !(CSobject@paired_end),
                                                          BPPARAM = bpParams)
        tsscounts.df <- SummarizedExperiment::assay(tsscounts)

        # sum all the TSS counts to gene counts
        id <- rownames(tsscounts.df)
        tsscounts.gene <-
            lapply(split(as.data.frame(tsscounts.df), id), colSums)
        tsscounts.gene <- do.call(rbind, tsscounts.gene)
        rm(tsscounts.df)

        # write back
        if (!is.na(outfile)) {
            write.table(
                tsscounts.gene,
                file = outfile,
                sep = "\t",
                quote = FALSE
            )
        }
        return(tsscounts.gene)

        }
)
