#' Plot read statistics from the CapSet object
#'
#' @param CSobject The \code{\link{CapSet}} object
#' @param plotType The type of plot to make. Choose from "stack" or "dodge" for either a stacked barchart,
#'                  or a bar chart with "dodged" positions (analogous to ggplot)
#' @param plotValue What values to plot. Choose from "numbers" or "proportions". If "proportions"
#'                  is selected, the proportion of reads w.r.t total demultiplexed reads per sample
#'                  would be plotted
#' @param outFile Output file name. (filename extention would be used to determine type).
#'                If outfile not specified, the plot would be retured on the screen
#'
#' @return A ggplot object, or a file. Plot showing the number/proportion of reads in each category, per sample
#' @export
#'
#' @examples
#'
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' plot_readStats(cs, plotType = "dodge", plotValue = "numbers", outFile = "test_numbers.pdf")
#'

plot_readStats <- function(CSobject,
                           plotType = c("stack", "dodge"),
                           plotValue = c("numbers", "proportions"),
                           outFile = NULL) {
    ## evaluate expressions
    stopifnot(is(CSobject, "CapSet"))
    stopifnot(plotType %in% c("stack", "dodge"))
    stopifnot(plotValue %in% c("numbers", "proportions"))

    ## get info on how many columns present
    si <- sampleInfo(CSobject)

    sicols <- colnames(si)
    fields_toplot <-
        c("demult_reads", "num_mapped", "num_filtered", "num_intss")
    msg <- "Plotting following information :"
    fields <- sapply(fields_toplot, function(x)
        if (x %in% sicols) {
            return(x)
        })
    message(cat(msg, fields))
    ## prepare df
    si_stats <- data.frame(
        sample = si$samples,
        demutiplexed_reads = si$demult_reads,
        mapped_reads = si$num_mapped
    )
    ## fill additional cols if present
    if (!(is.null(si$num_filtered)))
        si_stats$duplicate_free_reads <- si$num_filtered
    if (!(is.null(si$num_intss)))
        si_stats$reads_within_TSS <- si$num_intss

    if (plotValue == "proportions") {
        si_stats[-1] <- si_stats[-1] / si_stats$demutiplexed_reads
        # for stacked chart it's important to plot the cumulative difference of the numbers
        if (plotType == "stack") {
            si_stats[-1] <- get_stackedNum(si_stats[-1])
        }
        y_label <- "Proportion of demultiplexed reads"
    } else {
        if (plotType == "stack") {
            si_stats[-1] <- get_stackedNum(si_stats[-1])
        }
        y_label <- "Number of reads"
    }

    si_stats <- reshape2::melt(si_stats, id.vars = "sample")
    # plot stacked barchart
    p <-
        ggplot(si_stats, aes_string("sample", "value", fill = "variable")) +
        geom_bar(stat = "identity",
                 position = plotType,
                 color = "black") +
        theme_light(base_size = 16)  +
        scale_fill_brewer(type = "seq", palette = "YlGnBu") +
        coord_flip() +
        labs(x = "Sample", y = y_label, fill = "Category")
    # return
    if (!(is.null(outFile))) {
        ggsave(outFile, plot = p, dpi = 300)
    } else {
        return(p)
    }

}

## get cumulative differences of values from a DF, to plot stacked barchart of numbers
get_stackedNum <- function(df) {
    num <- apply(df, 1, function(x) {
        n <- c(abs(diff(as.numeric(x))), x[length(x)])
        return(n)
    })
    return(as.data.frame(t(num)))
}


#' Compare the precision of TSS detection between multiple samples
#'
#' @description Plot precision of TSS detection from multiple samples (bed files) with respect to
#' a given reference annotation.
#'
#' @param reference Reference Transcrips/Genes as a \code{\link{GRanges}} object
#' @param detectedTSS Either a CapSet object with TSS information (after running \code{\link{detect_TSS}}
#'                    or a character vector with paths to the BED files containing detcted TSSs
#' @param distanceCutoff Maximum distance (in base pairs) from reference TSS to plot
#' @param outFile Output file name (filename extention would be used to determine type)
#'                If outfile not specified, the plot would be returned on the screen
#' @param sampleNames Labels for input samples (in the same order as the input bed files)
#' @rdname plot_TSSprecision
#' @return A ggplot object, or a file. Plot showing perent of TSS detected per sample with respect to
#'         their cumulative distance to TSS of the provided reference
#'
#' @export
#' @examples
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' # load a txdb object
#' library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#' seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
#' transcripts <- transcripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#'
#' # Plotting the precision using a pre computed set of TSS (.bed files) :
#'
#' tssfile <- system.file("extdata", "testTSS_merged.bed", package = "icetea")
#' plot_TSSprecision(reference = transcripts, detectedTSS = tssfile,
#' 		sampleNames = "testTSS", distanceCutoff = 500,
#' 		outFile = "TSS_detection_precision.png")
#'

setMethod(
    plot_TSSprecision,
    signature = signature("GRanges", "character"),
    definition = function(reference,
                          detectedTSS,
                          distanceCutoff = 500,
                          outFile = NULL,
                          sampleNames) {
        # read bed files
        tssData <- lapply(detectedTSS, rtracklayer::import.bed)
        names(tssData) <- sampleNames
        # get plot
        plt <-
            plotPrecision(ref = reference,
                          tssData = tssData,
                          distCut = distanceCutoff)
        # output
        if (!(is.null(outFile))) {
            ggsave(outFile, plot = plt, dpi = 300)
        } else {
            return(plt)
        }
    }
)

#' Compare the precision of TSS detection between multiple samples
#'
#' @description Plot precision of TSS detection from multiple samples present within a
#' \code{\link{CapSet}} object, with respect to a given reference annotation.
#'
#' @docType methods
#' @rdname plot_TSSprecision
#' @param ... Additional arguments
#'
#' @export
#' @examples
#' ## Plotting the precision using a CapSet object :
#'
#' library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#' seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
#' # only use chrX to make the analysis faster
#' seqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "X"
#' transcripts <- transcripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#'
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' # plot
#' plot_TSSprecision(reference = transcripts, detectedTSS = cs,
#'                   outFile = "TSS_detection_precision.png")
#'

setMethod(
    plot_TSSprecision,
    signature = signature("GRanges", "CapSet"),
    definition = function(reference,
                          detectedTSS,
                          distanceCutoff = 500,
                          outFile = NULL,
                          ...) {
        # get the data out
        tssData <- detectedTSS@tss_detected

        if (is.null(tssData)) {
            stop("CapSet object does not contain the detected TSS information")
        }
        # get plot
        plt <-
            plotPrecision(ref = reference,
                          tssData = tssData,
                          distCut = distanceCutoff)
        # output
        if (!(is.null(outFile))) {
            ggsave(outFile, plot = plt, dpi = 300)
        } else {
            return(plt)
        }


    }
)

#' Plotprecision background script
#'
#' @param ref reference GRanges
#' @param tssData GRangesList object with TSS detected per sample
#' @param distCut max distance cutoff
#'
#' @importFrom ggplot2 aes_string stat_ecdf theme_light scale_x_continuous scale_color_brewer ggsave coord_flip
#' @return ggplot object
#'

plotPrecision <- function(ref, tssData, distCut) {
    # resize gene/transcript file to start
    refRanges <-
        GenomicRanges::resize(ref, width = 1, fix = "start")
    refRanges <- unique(refRanges)

    # get distances of bed entries from nearest TSS
    tssdistances <- lapply(tssData, function(x) {
        y <- GenomicRanges::distanceToNearest(x, refRanges)
        return(as.data.frame(y)$distance)
    })

    # melt to df
    tssdistances <- plyr::ldply(tssdistances, data.frame)
    colnames(tssdistances) <- c("sample", "distances")

    # plot ECDF with distance cutoff
    p <-
        ggplot(tssdistances, aes_string("distances", col = "sample")) +
        stat_ecdf(geom = "step", size = 1) +
        theme_light(base_size = 14)  +
        scale_x_continuous(limits = c(0, distCut)) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        scale_color_brewer(palette = "Set2") +
        labs(
            x = "Distances from nearby TSS (in bp)",
            y = "Cumulative Fraction",
            title = "TSS precisions",
            col = "Category"
        )
    return(p)

}
