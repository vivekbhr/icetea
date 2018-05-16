
#' Annotate the provided Transcription Start Sites
#'
#' This function annotates the provided TSS bed file to provide the number of TSS
#' falling within the genomic features from a given TxDB object. In order to break
#' ties between overlapping features, the function ranks the features by preference.
#' By default, the following order is used: fiveUTR > promoter > intron > coding >
#' spliceSite > threeUTR > intergenic. A custom order of feature ranks can also be
#' provided.
#'
#' @param tssBED A bed file with detected TSS/differential TSS coordinates
#' @param txdb A txdb object.
#' @param plotValue What values to plot (choose from "number", "percent" or NULL for no plot)
#' @param featureRank A vector with features to use for breaking ties, in decending order of
#'                    preference (highest to lowest),
#' @param outFile Output file name. (filename extention would be used to determine type).
#'                If outfile not specified, the plot would be retured on the screen
#'
#' @return A data.frame with number of TSS falling into each feature
#' @import TxDb.Dmelanogaster.UCSC.dm6.ensGene
#' @importFrom ggplot2 ggplot aes_string geom_bar scale_fill_brewer labs theme
#'                     theme_gray coord_flip ggsave
#' @importFrom stats reshape
#' @importFrom VariantAnnotation locateVariants AllVariants PromoterVariants
#' @export
#'
#' @examples
#' # load a txdb object
#' library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#' seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
#' # limiting the annotation to X chromosome
#' seqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "X"
#'
#' # annotate a given TSS bed file
#' dir <- system.file("extdata", package = "icetea")
#' tss <- file.path(dir, "testTSS_merged.bed")
#' annotations <- annotateTSS(tssBED = tss, TxDb.Dmelanogaster.UCSC.dm6.ensGene,
#'                plotValue = "number", outFile = "TSS_annot.pdf")
#'

annotateTSS <- function(tssBED,
                        txdb,
                        featureRank = c("fiveUTR",
                                        "promoter",
                                        "intron",
                                        "coding",
                                        "spliceSite",
                                        "threeUTR",
                                        "intergenic"),
                         plotValue = "number",
                         outFile = NULL) {
    ## resolve 1:many mapping issue by prioritising some features over others
    stopifnot(length(featureRank) == 7)
    rankvec <- seq_len(7)
    names(rankvec) <- featureRank

    # get data
    tss <- rtracklayer::import.bed(tssBED)
    # Annotate
    suppressWarnings({
        db <- locateVariants(
        query = tss,
        subject = txdb,
        AllVariants(promoter = PromoterVariants(
            upstream = 500, downstream = 0
            ))
        )
    })
    ## resolve 1:many mapping isues using ranks from rankdf
    df <- data.frame(QUERYID = db$QUERYID,
                    LOCATION = db$LOCATION)
    df$rank <- vapply(as.character(df$LOCATION), getranks, rank_vec = rankvec, FUN.VALUE = numeric(length = 1))
    df2 <- splitranks(df)
    ## Return a table of tss counts per feature
    final_table <- as.data.frame(table(df2$LOCATION))
    colnames(final_table) <- c("feature", "value")
    ## plot if asked
    if (!is.null(plotValue)) {
        if (plotValue == "number") {
            n <- "Number "
        } else if (plotValue == "percent") {
            final_table$value <- (final_table$value / sum(final_table$value)) * 100
            n <- "% "
        } else {
            warning("Plot type neither 'number' nor 'percent'.")
        }

        p <- ggplot(final_table,
               aes_string("feature", "value", fill = "feature")) +
            geom_bar(stat = "identity", position = "dodge") +
            scale_fill_brewer(palette = "Set1") +
            labs(x = "Feature", y = paste0(n, "of TSS")) +
            theme(legend.position = "none") +
            theme_gray(base_size = 16) +
            coord_flip()

        if (is.null(outFile)) {
            print(p)
        } else {
            ggsave(plot = p, filename = outFile)
        }
    }

    return(final_table)

}


#' Assign feature ranks on a VariantAnnotation output
#'
#' @param x output from VariantAnnotation
#' @param rank_vec the pre-set vector of ranks
#'
#' @return A vector of ranks of length = length of input features
#'
#'
getranks <- function(x, rank_vec) {
                    return(rank_vec[names(rank_vec) == x])
                }

#' Get features with the best rank for each TSS
#'
#' @param x output of getranks
#'
#' @return A data frame with counts
#'
#'
splitranks <- function(x) {
    l <- lapply(split(x, x$QUERYID), unique)
    l2 <- lapply(l, function(y) {
        return(y[which(y$rank == min(y$rank)), ])
    })
    l3 <- do.call(rbind, l2)
    return(l3)
}
