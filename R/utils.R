#' Count the number of reads in a given GRanges
#'
#' @param regions The GRanges object
#' @param bams path to bam files from where the reads have to be counted
#'
#' @return Total counts within given ranges per BAM file.
#'
numReadsInBed <- function(regions, bams = NA) {
    counts <-
        GenomicAlignments::summarizeOverlaps(
            GenomicRanges::GRangesList(regions),
            reads = Rsamtools::BamFileList(as.character(bams)),
            mode = "Union",
            inter.feature = FALSE
        )
    numreads <- SummarizedExperiment::assay(counts)
    return(t(numreads))
}


#' Match BAM headers bw files and get active chromosome list (from restrict)
#' (written by Aaron Lun, 12 Dec 2014, copied and modified here)
#'
#' @param bam.files Character vector (bam files)
#' @param restrict Chromosomes to select
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
#' @param bamFiles Character vector (bam files)
#' @param restrictChr Chromosomes to select
#' @param binSize Size of bins
#'
#' @return GRanges (bins) for both strands
#'
getChromBins <- function(bamFiles, restrictChr = NULL, binSize) {

    keptChrs <- activeChrs(bamFiles, restrict = restrictChr)
    instances <- as.integer(keptChrs/binSize)

    gr.bins.plus <- lapply(names(keptChrs), function(x){
        gr <- GenomicRanges::GRanges(seqnames = x,
                        ranges = IRanges::successiveIRanges(rep(binSize, keptChrs[x])),
                        strand = "+"
                        )
        return(gr)
        })
    gr.bins.plus <- suppressWarnings({do.call("c", gr.bins.plus)})
    GenomeInfoDb::seqlengths(gr.bins.plus) <- keptChrs
    gr.bins.minus <- gr.bins.plus
    GenomicRanges::strand(gr.bins.minus) <- "-"
    return(list(gr.plus = GenomicRanges::trim(gr.bins.plus),
                gr.minus = GenomicRanges::trim(gr.bins.minus))
    )
}

#' preprocess reads to count only 5' overlaps
#'
#' @param reads GAlignment object
#' @param width New read length
#' @param fix 'Start' for 5'
#' @param ... Other
#'
#' @return Resized reads as GRanges
#'
ResizeReads <- function(reads, width = 1, fix = "start", ...) {
    reads <- as(reads, "GRanges")
    stopifnot(all(GenomicRanges::strand(reads) != "*"))
    GenomicRanges::resize(reads, width = width, fix = fix, ...)
}

#' Calculate local enrichment of windows over background
#' 'local' filter copied and modified from csaw::filterWindows
#' written by Aaron Lun (5 November 2014, last modified 3 March 2017)
#'
#' @param data RangedSE object (windows)
#' @param background RangedSE object (background)
#' @param assay.data Arg to pass forward
#' @param assay.back Arg to pass forward
#'
#' @return list with data and background abundances
#'
localFilter <- function(data,
                        background,
                        assay.data = 1,
                        assay.back = 1) {
    # get default prior counts
    prior.count <- formals(edgeR::aveLogCPM.DGEList)$prior.count

    # no need to get "effective width", fraglength is always 1 for 5' counting
    dwidth <- GenomicRanges::width(SummarizedExperiment::rowRanges(data))
    bwidth <- GenomicRanges::width(SummarizedExperiment::rowRanges(background))
    # avg logCPM of data
    abundances <- csaw::scaledAverage(
                        csaw::asDGEList(
                            data, assay = assay.data),
                        scale = 1, prior.count = prior.count)

    if (!identical(nrow(data), nrow(background))) {
        stop("data and background should be of the same length")
    }

    relative.width <- (bwidth  - dwidth)/dwidth
    bg.y <- csaw::asDGEList(background, assay = assay.back)
    bg.y$counts <- bg.y$counts - SummarizedExperiment::assay(data, assay = assay.data)

    # Some protection for negative widths
    # (counts should be zero, so only the prior gets involved in bg.ab).
    subzero <- relative.width <= 0
    if (any(subzero)) {
        relative.width[subzero] <- 1
        bg.y$counts[subzero,] <- 0L
    }
    # avg logCPM of background
    bg.ab <- csaw::scaledAverage(bg.y, scale = relative.width, prior.count = prior.count)
    # filter stat is the fold change of data over bg
    filter.stat <- abundances - bg.ab

    return(list(abundances = abundances,
                back.abundances = bg.ab,
                filter = filter.stat)
           )
}
