#' Get gene-level counts from TSS data
#'
#' @param CSobject The \code{\link{CapSet}} object to use.
#' @param transcriptGRL A GRangesList object containing transcripts, created using transcriptsBy(txdb)
#' @param regionAroundTSS How many bases downstream of TSS to count
#' @param single_end Logical, indicating whether reads are single end
#' @param outfile Tab-separated output file name (if required)
#'
#' @return data.frame with gene-level counts for all genes in the txdb object
#'
#' @importFrom GenomicFeatures transcriptsBy
#' @export
#'
#' @examples
#'
#'  # load a txdb object
#'  library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#'  seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
#'  # get transcripts by gene
#'  dm6trans <- transcriptsBy(TxDb.Dmelanogaster.UCSC.dm6.ensGene, "gene")
#'  # load a CapSet object
#'  cs <- exampleCSobject()
#'  # get gene counts, counting reads around 500 bp of the TSS
#'  gcounts <- get_geneCounts(cs, dm6trans)
#'

get_geneCounts <- function(CSobject, transcriptGRL, regionAroundTSS = 500, single_end = TRUE, outfile = NA) {

    # add gene names and unlist
    transcriptGRL <- mapply(function(x, name) {
                                    x$geneID <- name
                                    return(x)
                                    }, transcriptGRL, names(transcriptGRL))

    # get XX bp region around transcripts to count the reads
    transcriptGR <- unlist(GenomicRanges::GRangesList(transcriptGRL))
    transcriptGR <- GenomicRanges::resize(transcriptGR, regionAroundTSS, fix = "start")
    # get bamfiles
    si <- sampleInfo(CSobject)
    bamfiles <- si$filtered_file

    # count reads
    tsscounts <- GenomicAlignments::summarizeOverlaps(transcriptGR,
                                                reads = Rsamtools::BamFileList(bamfiles),
                                                singleEnd = single_end)
    tsscounts.df <- SummarizedExperiment::assay(tsscounts)

    # sum all the TSS counts to gene counts
    id <- rownames(tsscounts.df)
    tsscounts.gene <- lapply(split(as.data.frame(tsscounts.df), id), colSums)
    tsscounts.gene <- do.call(rbind, tsscounts.gene)
    rm(tsscounts.df)

    # write back
    if(!is.na(outfile)) {
    write.table(tsscounts.gene,
        file = outfile,
        sep = "\t", quote = FALSE)
    }
    return(tsscounts.gene)

}
