#' Get platform-specific multicore params
#'
#' @param cores integer. No. of cores to use.
#' @importFrom BiocParallel SnowParam MulticoreParam
#'
#' @return BPPARAM object
#'
getMCparams <- function(cores) {
    if (cores == 1) {
        param <- BiocParallel::SerialParam()
    } else {
        param <- switch(Sys.info()[['sysname']],
                        Windows = {return(SnowParam(workers = cores))},
                        Linux = {return(MulticoreParam(workers = cores))},
                        Darwin = {return(MulticoreParam(workers = cores))}
                    )
    }
    return(param)
}

#' Get flags to read from bam
#'
#' @param countAll logical. count all reads?
#' @importFrom Rsamtools scanBamFlag
#'
#' @return bamFlags
#'
getBamFlags <- function(countAll) {
        # get baminfo
    if (isTRUE(countAll)) {
        # if countAll given, count both reads (in PE mode) or all reads (in SE mode)
        bamFlags <- scanBamFlag(
                                isUnmappedQuery = FALSE,
                                isSecondaryAlignment = FALSE
                            )
    } else {
        # else count only R1 (in PE mode)
        bamFlags <- scanBamFlag(
                                isUnmappedQuery = FALSE,
                                isFirstMateRead = TRUE,
                                isSecondaryAlignment = FALSE
                            )
    }
    return(bamFlags)
}

#' Count the number of reads in a given GRanges
#'
#' @param regions The GRanges object to count reads in.
#' @param bams character. path to bam files from where the reads have to be counted
#' @param countall logical. whether to keep both reads of paired-end data
#'
#' @return Total counts within given ranges per BAM file.
#'
numReadsInBed <- function(regions, bams = NA, countall = FALSE) {
    counts <-
        GenomicAlignments::summarizeOverlaps(
            GenomicRanges::GRangesList(regions),
            reads = Rsamtools::BamFileList(as.character(bams)),
            mode = "Union",
            inter.feature = FALSE,
            param = Rsamtools::ScanBamParam(flag = getBamFlags(countAll = countall))
        )
    numreads <- SummarizedExperiment::assay(counts)
    return(t(numreads))
}

#' Match BAM headers bw files and get active chromosome list (from restrict)
#' (written by Aaron Lun, 12 Dec 2014, copied and modified here)
#'
#' @param bam.files Character . bam files to check
#' @param restrict character. Chromosomes to select
#'
#' @return Vector of selected chromosomes
#'
activeChrs <- function(bam.files, restrict)
{
    keptChrs <- NULL
    for (bam in bam.files) {
        chrs <- Rsamtools::scanBamHeader(bam)[[1]][[1]]
        chrs <- chrs[order(names(chrs))]
        if (is.null(keptChrs)) { keptChrs <- chrs }
        else if (!identical(keptChrs, chrs)) {
            warning("chromosomes are not identical between BAM files")
            pairing <- match(names(keptChrs), names(chrs))
            keptChrs <- pmin(keptChrs[!is.na(pairing)], chrs[pairing[!is.na(pairing)]])
        }
    }
    if (length(restrict)) { keptChrs <- keptChrs[names(keptChrs) %in% restrict] }
    return(keptChrs)
}

#' Get chromosome bins from BAM files
#'
#' @param bamFiles Character. bam files
#' @param restrictChr character. Chromosomes to select
#' @param binSize numeric. Size of bins
#'
#' @importFrom GenomicRanges tileGenome
#' @return GRanges (bins) for both strands
#'
getChromBins <- function(bamFiles, restrictChr = NULL, binSize) {
    keptChrs <- activeChrs(bamFiles, restrict = restrictChr)
    gr.bins.plus <- tileGenome(keptChrs, tilewidth = 10, cut.last.tile.in.chrom = TRUE)
    gr.bins.minus <- gr.bins.plus
    GenomicRanges::strand(gr.bins.plus) <- "+"
    GenomicRanges::strand(gr.bins.minus) <- "-"
    return(list(gr.plus = gr.bins.plus,
                gr.minus = gr.bins.minus))
}

#' Get chromosome sliding windows from BAM files
#'
#' @param bamFiles Character vector (bam files)
#' @param restrictChr Chromosomes to select
#' @param binSize Size of bins
#' @param stepSize Size of window slide
#'
#' @return GRanges (sliding windows) for both strands
#'
getChromWindows <- function(bamFiles, restrictChr = NULL, binSize, stepSize) {
    keptChrs <- activeChrs(bamFiles, restrict = restrictChr)
    gr.total <- GenomicRanges::GRanges(
                                        seqnames = names(keptChrs),
                                        ranges = IRanges::IRanges(start = 1, end = keptChrs),
                                        strand = "+")
    gr.bins.plus <- GenomicRanges::slidingWindows(gr.total, width = binSize, step = stepSize)
    gr.bins.plus <- unlist(gr.bins.plus)
    gr.bins.minus <- gr.bins.plus
    GenomicRanges::strand(gr.bins.minus) <- "-"
    return(list(gr.plus = gr.bins.plus,
                gr.minus = gr.bins.minus))
}


#' preprocess reads to count only 5' overlaps
#'
#' @param reads GAlignment object to resize
#' @param width integer. New read length
#' @param fix character. 'Start' for 5'
#'
#' @return Resized reads as GRanges
#'
readsTo5p <- function(reads, width = 1, fix = "start") {
    reads <- as(reads, "GRanges")
    stopifnot(all(GenomicRanges::strand(reads) != "*"))
    GenomicRanges::resize(reads, width = width, fix = fix)
}

#' preprocess reads to count only 3' overlaps
#'
#' @param reads GAlignment object to resize
#' @param width integer. New read length
#' @param fix character. 'Start' for 5'
#'
#' @return Resized reads as GRanges
#'
readsTo3p <- function(reads, width = 1, fix = "end") {
    reads <- as(reads, "GRanges")
    stopifnot(all(GenomicRanges::strand(reads) != "*"))
    GenomicRanges::resize(reads, width = width, fix = fix)
}

#' preprocess reads to count only center overlaps
#'
#' @param reads GAlignment object to resize
#' @param width integer. New read length
#' @param fix character. 'Start' for 5'
#'
#' @return Resized reads as GRanges
#'
readsToCenter <- function(reads, width = 1, fix = "center") {
    reads <- as(reads, "GRanges")
    stopifnot(all(GenomicRanges::strand(reads) != "*"))
    GenomicRanges::resize(reads, width = width, fix = fix)
}


# Calculate local enrichment of windows over background
# 'local' filter copied and modified from csaw::filterWindows
# written by Aaron Lun (5 November 2014, last modified 3 March 2017)
#
# @param data RangedSE object (windows)
# @param background RangedSE object (background)
# @param assay.data Arg to pass forward
# @param assay.back Arg to pass forward
#
# @return list with data and background abundances
#
localFilter <- function(data,
                        background,
                        assay.data = 1,
                        assay.back = 1) {
    # check lengths
    if (!identical(nrow(data), nrow(background))) {
        stop("data and background should be of the same length")
    }
    # get default prior counts
    prior.count <- formals(edgeR::aveLogCPM.DGEList)$prior.count

    # no need to get "effective width", fraglength is always 1 for 5' counting
    dwidth <- GenomicRanges::width(SummarizedExperiment::rowRanges(data))
    bwidth <- GenomicRanges::width(SummarizedExperiment::rowRanges(background))

    # avg logCPM of data
    #    dat.y <- csaw::asDGEList(data, assay = assay.data)
    #    dat.y <- edgeR::estimateCommonDisp(dat.y)
    data.ab <- csaw::scaledAverage(data, scale = 1,
                                   assay.id = 1,
                                   prior.count = prior.count)

    relative.width <- (bwidth  - dwidth)/dwidth
    #    bg.y <- csaw::asDGEList(background, assay = assay.back)
    #    bg.y <- edgeR::estimateCommonDisp(bg.y)
    #    bg.y$counts <- bg.y$counts - dat.y$counts
    bg.y <- background
    SummarizedExperiment::assay(bg.y) <-
                   SummarizedExperiment::assay(bg.y) -
                   SummarizedExperiment::assay(data)

    # Some protection for negative widths
    # (counts should be zero, so only the prior gets involved in bg.ab).
    subzero <- relative.width <= 0
    if (any(subzero)) {
        relative.width[subzero] <- 1
        bg.y$counts[subzero,] <- 0L
    }
    # avg logCPM of background
    bg.ab <- csaw::scaledAverage(bg.y, scale = relative.width,
                                 assay.id = 1,
                                 prior.count = prior.count)
    # filter stat is the fold change of data over bg
    filter.stat <- data.ab - bg.ab

    return(list(abundances = data.ab,
                back.abundances = bg.ab,
                logFC = filter.stat)
    )
}
