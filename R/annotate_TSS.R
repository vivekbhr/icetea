
#' Annotate detected Transcription Start Sites
#'
#' @param tssGR GRanges object with detected TSS/differential TSS results.
#' @param txdb A txdb object.
#' @param plot Type of plot to make (choose from "number", "percent" or NA for no plot)
#'
#' @return Annotation of detected TSS
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_bar scale_fill_brewer labs theme theme_gray
#' @importFrom stats reshape
#'
#' @examples
#' # load a txdb object
#' library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#' seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
#' # annotate a given TSS bed file
#' tss <- load(file.path(dir, "testTSS_grl.Rdata"))
#'
#' annotations <- annotate_TSS(tssGR = tss$wt, TxDb.Dmelanogaster.UCSC.dm6.ensGene, plot = "testplot.pdf")
#'

annotate_TSS <- function(tssGR, txdb, plot = NA) {

    ## resolve 1:many mapping issue by prioritising some features over others
    rank_df <- data.frame(feature = c("fiveUTR", "promoter", "intron", "coding",
      "spliceSite", "threeUTR", "intergenic"),
          rank = c(1,2,3,4,5,6,7))
    # Annotate
    db <- VariantAnnotation::locateVariants(query = tssbed,
                                            subject = txdb,
    VariantAnnotation::AllVariants(
    promoter = VariantAnnotation::PromoterVariants(upstream = 500, downstream = 0)))
    ## resolve 1:many mapping isues using ranks from rankdf
    t <- data.frame(QUERYID = db$QUERYID, LOCATION = db$LOCATION)
    tt <- getranks(t, rankdf = rank_df)
    ttt <- splitranks(tt)
    ## Return a table of tss counts per feature
    final_table <- as.data.frame(table(ttt$LOCATION))
    colnames(final_table) <- c("feature", "value")
    ## plot if asked
    if (!is.na(plot)) {
    if (plot == "number") {
    n <- "Number "
    } else if(plot == "percent") {
    final_table$value <- (final_table$value/sum(final_table$value))*100
    n <- "% "
    } else {
    warning("Plot type neither 'number' nor 'percent'.")
    }

    print(ggplot(final_table, aes_string("feature", "value", fill = "feature")) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = "Set1") +
    labs(x = "Feature", y = paste0(n, "of TSS")) +
    theme(legend.position = "none") +
    theme_gray(base_size = 16)
    )
    }

    return(final_table)

}


#' Assign feature ranks on a VariantAnnotation output
#'
#' @param x output from VariantAnnotation
#' @param rankdf the defined rankdf data.frame
#'
#' @return A list of ranks
#'
#'
getranks <- function(x, rankdf) {
    x$rank <- sapply(x$LOCATION, function(y) {
    return(rankdf[rankdf$feature == y, "rank"])
    })
    return(x)
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
    return(y[which(y$rank == min(y$rank)),])
    })
    l3 <- plyr::ldply(l2, data.frame)
    return(l3)
}


# Melt the output df from splitranks
melt <- function(x) {
    vars <- colnames(x[2:ncol(x)])
    d <- plyr::unrowname(reshape(x, direction = "long", idvar = ".id", varying = vars) )
    d$time <- vars[d$time]
    colnames(d) <- c("variable", "Feature", "value")
}
