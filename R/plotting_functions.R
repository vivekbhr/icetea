#' Plot read statistics from the CapSet object
#'
#' @rdname plotReadStats
#' @param CSobject The \code{\link{CapSet}} object
#' @param plotType character. The type of plot to make. Choose from "stack" or "dodge" for either a stacked barchart,
#'                  or a bar chart with "dodged" positions (analogous to ggplot)
#' @param plotValue character. What values to plot. Choose from "numbers" or "proportions". If "proportions"
#'                  is selected, the proportion of reads w.r.t total demultiplexed reads per sample
#'                  would be plotted
#' @param outFile character. Output file name. (filename extention would be used to determine type).
#'                If outfile not specified, the plot would be retured on the screen
#'
#' @return A ggplot object, or a file. Plot showing the number/proportion of reads in each category, per sample
#'
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_light scale_fill_brewer coord_flip labs ggsave
#' @export
#'
#' @examples
#'
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' plotReadStats(cs, plotType = "dodge", plotValue = "numbers", outFile = "test_numbers.pdf")
#'
setMethod(
    "plotReadStats",
    signature = "CapSet",
    definition = function(
                        CSobject,
                        plotType,
                        plotValue,
                        outFile) {
    ## evaluate expressions
    stopifnot(plotType %in% c("stack", "dodge"))
    stopifnot(plotValue %in% c("numbers", "proportions"))
    fields_toplot <- c("demult_reads", "num_mapped", "num_filtered", "num_intss")

    ## get info on how many columns to plot
    si <- sampleInfo(CSobject)[c("samples",fields_toplot)]
    # subset si for non-na cols
    nacols <- apply(si, 2, function(x) all(is.na(x)))
    si <- si[!nacols]
    msg <- "Plotting following information :"
    message(cat(msg, colnames(si) ))

    ## prepare df
    si_stats <- data.frame(sample = si$samples)

    ## fill additional cols if present
    if (!is.null(si$demult_reads)) si_stats$demultiplexed_reads <- si$demult_reads
    if (!is.null(si$num_mapped)) si_stats$mapped_reads <- si$num_mapped
    if (!is.null(si$num_filtered)) si_stats$duplicate_free_reads <- si$num_filtered
    if (!is.null(si$num_intss)) si_stats$reads_within_TSS <- si$num_intss

    if (plotValue == "proportions") {
        # plot proportion w.r.t lowst category
        basecat <- "demultiplexed_reads"
        if (is.null(si_stats[, basecat])) basecat <- "mapped_reads"
        if (is.null(si_stats[, basecat])) basecat <- "duplicate_free_reads"
        if (is.null(si_stats[, basecat])) stop("Can't plot proportions with only one category")

        si_stats[-1] <- si_stats[-1] / si_stats[, basecat]
        # for stacked chart it's important to plot the cumulative difference of the numbers
        if (plotType == "stack") {
            si_stats[-1] <- get_stackedNum(si_stats[-1])
        }
        y_label <- paste0("Proportion of ", basecat, " reads")

    } else if (plotValue == "numbers") {
        if (plotType == "stack") {
            si_stats[-1] <- get_stackedNum(si_stats[-1])
        }
        y_label <- "Number of reads"
    }

    varcols <- colnames(si_stats)[2:ncol(si_stats)]
    si_stats <- reshape(si_stats, direction = "long", idvar = "sample",
                        varying = varcols, timevar = "variable",
                        times = varcols, v.names = "value")
    rownames(si_stats) <- NULL
    si_stats$variable <- factor(si_stats$variable,
                                levels = c("demultiplexed_reads",
                                           "mapped_reads",
                                           "duplicate_free_reads",
                                           "reads_within_TSS"))
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

})

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
#' @param detectedTSS Either a CapSet object with TSS information (after running \code{\link{detectTSS}}
#'                    or a character vector with paths to the BED files containing detcted TSSs
#' @param distanceCutoff integer. Maximum distance (in base pairs) from reference TSS to plot
#' @param outFile character. Output file name (filename extention would be used to determine type)
#'                If outfile not specified, the plot would be returned on the screen
#' @param sampleNames character. Labels for input samples (in the same order as the input bed files)
#' @rdname plotTSSprecision
#' @return A ggplot object, or a file. Plot showing perent of TSS detected per sample with respect to
#'         their cumulative distance to TSS of the provided reference
#'
#' @export
#' @examples
#'
#' # load a txdb object
#' library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#' seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
#' transcripts <- transcripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#'
#' # Plotting the precision using a pre computed set of TSS (.bed files) :
#'
#' tssfile <- system.file("extdata", "testTSS_merged.bed", package = "icetea")
#' plotTSSprecision(reference = transcripts, detectedTSS = tssfile,
#'                 sampleNames = "testTSS", distanceCutoff = 500,
#'                 outFile = "TSS_detection_precision.png")
#'

setMethod(
    plotTSSprecision,
    signature = signature("GRanges", "character"),
    definition = function(
                        reference,
                        detectedTSS,
                        distanceCutoff = 500,
                        outFile = NULL,
                        sampleNames) {
        # read bed files as GRangesList
        tssData <- GenomicRanges::GRangesList(
            lapply(detectedTSS, rtracklayer::import.bed) )
        names(tssData) <- sampleNames
        # get plot
        plt <- plotPrecision(ref = reference,
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
#' @rdname plotTSSprecision
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
#' plotTSSprecision(reference = transcripts, detectedTSS = cs,
#'                   outFile = "TSS_detection_precision.png")
#'

setMethod(
    plotTSSprecision,
    signature = signature("GRanges", "CapSet"),
    definition = function(reference,
                          detectedTSS,
                          distanceCutoff = 500,
                          outFile = NULL,
                          ...) {
        # get the data out
        tssData <- GenomicRanges::GRangesList(detectedTSS@tss_detected)

        if (is.null(tssData)) {
            stop("CapSet object does not contain the detected TSS information")
        }
        # get plot
        plt <- plotPrecision(
                            ref = reference,
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
#' @param ref GRanges. reference ranges to compare the precision with.
#' @param tssData GRangesList object with TSS detected per sample
#' @param distCut integer. max distance cutoff
#'
#' @importFrom ggplot2 aes_string stat_ecdf theme_light scale_x_continuous
#'             scale_y_continuous scale_color_brewer ggsave coord_flip
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
    samplenames <- vapply(tssdistances, length, integer(1L))
    samplenames <- rep(names(samplenames), samplenames)
    sampledists <- do.call(c, tssdistances)
    tssdistances <- data.frame(sample = samplenames, distances = sampledists)
    #colnames(tssdistances) <- c("sample", "distances")
    # print message for removed values
    removed <- tssdistances$distances > distCut
    tssdistances[removed,] -> highdist
    vapply(split(highdist, factor(highdist$sample)), nrow, numeric(1L)) -> removedNums

    message(paste0("There are ", sum(removed),
                    " regions with distance > ", distCut,
                    " bp to the closest TSS. They are all being reduced to ", distCut,
                    " bp for the calculation. Samplewise numbers are : ", removedNums))
    tssdistances$distances[removed] <- distCut
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
            col = "Sample"
        )
    return(p)

}
