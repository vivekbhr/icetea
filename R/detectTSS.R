#' Perform stranded Bin counts
#'
#' @param bam.files character vector. BAM files to use
#' @param restrictChrs character vector. chromosomes to use
#' @param bam_param ScanBAMParams
#' @param bp_param BPPARAM
#' @param window_size integer. size of window to use
#' @param sliding logical. perform sliding window counts?
#'
#' @importFrom SummarizedExperiment assay rowRanges
#'
#' @return RangedSE object with forward and reverse strand counts
#'
strandBinCounts <- function(bam.files, restrictChrs, bam_param, bp_param, window_size, sliding = FALSE) {

    if (sliding == FALSE) {
        windows <- getChromBins(bam.files, restrictChr = restrictChrs, binSize = window_size)
        ignoreMultiMap <- TRUE
    } else {
        windows <- getChromWindows(bam.files, restrictChr = restrictChrs,
                                    binSize = window_size, stepSize = floor(window_size/2) )
        ignoreMultiMap <- FALSE
    }

    fdata <-
        GenomicAlignments::summarizeOverlaps(
            features = windows$gr.plus,
            reads = bam.files,
            mode = "IntersectionStrict",
            ignore.strand = FALSE,
            inter.feature = ignoreMultiMap,
            singleEnd = TRUE,
            fragments = FALSE,
            preprocess.reads = ResizeReads,
            param = bam_param,
            BPPARAM = bp_param)

    rdata <-
        GenomicAlignments::summarizeOverlaps(
            features = windows$gr.minus,
            reads = bam.files,
            mode = "IntersectionStrict",
            ignore.strand = FALSE,
            inter.feature = ignoreMultiMap,
            singleEnd = TRUE,
            fragments = FALSE,
            preprocess.reads = ResizeReads,
            param = bam_param,
            BPPARAM = bp_param)
    coldat <- S4Vectors::DataFrame(bam.files = bam.files,
                                    forward.totals = S4Vectors::colSums(assay(fdata)),
                                    reverse.totals = S4Vectors::colSums(assay(rdata)),
                                    ext = NA,
                                    rlen = 1L)
    combined <- SummarizedExperiment::SummarizedExperiment(
                            rbind(assay(fdata, "counts"), assay(rdata, "counts")),
                            rowRanges = c(rowRanges(fdata), rowRanges(rdata)),
                            colData = coldat)
    # drop empty bins
    combined <- combined[BiocGenerics::rowSums(assay(combined)) > 0]
    # Suggestion : Drop bins with counts < threshold ?
    combined$totals <- combined$forward.totals + combined$reverse.totals
    return(combined)
}

#' Detection of Trancription start sites based on local enrichment
#'
#' @rdname detectTSS
#' @param CSobject CapSet object created using \code{\link{newCapSet}} function
#' @param groups a character vector that contains group name of the sample, for replicate-based TSS
#'               calling (see example)
#' @param outfile_prefix Output name prefix for the .Rdata file containing window counts, background counts
#'                       and filtering statistics calculated during TSS detection.
#' @param windowSize Size of the window to bin the genome for TSS detection. By default, a window size of
#'                   10 is used for binning the genome, however smaller window sizes can optionally be provided
#'                   for higher resolution TSS detection. Note that the background size is set to 200x the
#'                   window size (2kb for 10bp windows) to calculate local enrichment. Adjacent enriched windows
#'                   are merged with a distance cutoff, which is the same as window size to get final TSS widths.
#' @param sliding TRUE/FALSE. Indicating whether or not to use sliding windows. The windows are shifted by length which
#'                is half of the specified window length.
#' @param foldChange Numeric. A fold change cutoff of local enrichment to detect the TSS. For samples with
#'        usual' amount of starting material and squencing depth (>=5ug starting material,
#'        = 5 mil reads/sample), a cut-off of 6 fold can be used. For samples with low
#'        amount of material or sequencing depth, use a lower cut-off (eg. use 2-fold for
#'        samples with 500ng starting material).
#' @param mergeLength Integer. Merge the windows within this distance that pass the foldChange cutoff.
#'                    Default (1L) means that only subsequently enriched windows would be merged.
#' @param restrictChr Chromosomes to restrict the analysis to.
#' @param ncores No. of cores/threads to use
#'
#' @return .bed files containing TSS position for each group, along with a bed file for consensus
#'        (union) TSS sites of all samples.
#'
#' @export
#' @importFrom utils write.table
#' @importFrom SummarizedExperiment mcols mcols<- colData colData<- rowRanges
#' @importFrom csaw readParam strandedCounts regionCounts filterWindows mergeWindows
#'
#' @examples
#'
#' # before running this
#' # 1. Create a CapSet object
#' # 2. de-multiplex the fastqs
#' # 3. map them
#' # 4. filter duplicate reads from mapped BAM
#'
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' # detect TSS (samples in same group are treated as replicates)
#' cs <- detectTSS(cs, groups = rep(c("wt","mut"), each = 2), outfile_prefix = "testTSS",
#'            foldChange = 6, restrictChr = "X", ncores = 1)
#'

setMethod("detectTSS",
          signature = "CapSet",
          function(CSobject,
                    groups,
                    outfile_prefix,
                    windowSize,
                    sliding,
                    foldChange,
                    mergeLength,
                    restrictChr,
                    ncores
                    ) {
            # check whether group and outfile_prefix is provided
            if (missing(outfile_prefix))
                stop("Please provide outfile_prefix!")
            if (missing(groups))
                stop("Please provide groups!")

            # convert group to char
            si <- sampleInfo(CSobject)
            design <-
                data.frame(row.names = si$samples, group = as.character(groups))

            if (all(is.na(si$filtered_file))) {
                warning("Filtered files not found under sampleInfo(CSobject). Using mapped files")
                bam.files <- si$mapped_file
            } else {
                bam.files <- si$filtered_file
            }
            if (any(is.na(bam.files))) stop("Some or all of the bam files are not defined!")
            if (sum(file.exists(bam.files)) != length(bam.files)) {
                stop("One or more bam files don't exist! Check sampleInfo(CSobject) ")
            }

            # Counting params
            countall = !(CSobject@paired_end)
            bamParams <- Rsamtools::ScanBamParam(
                                flag = getBamFlags(countAll = countall))

            bpParams <- getMCparams(ncores)
            # register parallel backend
            if (!BiocParallel::bpisup(bpParams)) {
                BiocParallel::bpstart(bpParams)
                on.exit(BiocParallel::bpstop(bpParams))
            }
            # window size
            bin_size <- windowSize
            # background region size (200 x Window size)
            surrounds <- 200*bin_size

            # Count reads into sliding windows
            data <- strandBinCounts(bam.files, restrictChr,
                                    bam_param = bamParams,
                                    bp_param = bpParams,
                                    window_size = bin_size,
                                    sliding = sliding)
              # add metadata
            #mdat <- list(spacing = bin_size, width = bin_size,
            #            shift = 0, bin = TRUE, final.ext = 1)
            #S4Vectors::metadata(data) <- mdat
            colnames(data) <- rownames(design)
            colData(data) <- c(colData(data), design)

            # Get counts for background region
            neighbors <- suppressWarnings(GenomicRanges::trim(
                GenomicRanges::resize(rowRanges(data),
                                        surrounds, fix = "center")
            ))

            wider <-
                suppressWarnings({
                GenomicAlignments::summarizeOverlaps(
                    features = neighbors,
                    reads = bam.files,
                    mode = "IntersectionStrict",
                    ignore.strand = FALSE,
                    inter.feature = FALSE,
                    singleEnd = TRUE,
                    fragments = FALSE,
                    preprocess.reads = ResizeReads,
                    param = bamParams,
                    BPPARAM = bpParams)
                  })

            #S4Vectors::metadata(wider) <- mdat
            # set totals to same value as data (to avoid error from filterWindows)
            colData(wider) <- colData(data)
            colnames(wider) <- rownames(design)
            colData(wider) <- c(colData(wider), design)

            ## take out groups --> Generate filter statistics for each group (based on local enrichment)
            filterstat <- lapply(unique(design$group), function(x) {
                stat <- localFilter(data[, data$group == x],
                                        wider[, wider$group == x])
                  return(S4Vectors::DataFrame(stat))
            })
            # add filter stats as metadata to the data
            mcols(data) <- filterstat

            # Require X-fold enrichment over local background to keep the window (similar to MACS)
            keep <- lapply(filterstat, function(x) {
                kp <- x$filter > log2(foldChange)
                return(kp)
            })

            filtered.data <- lapply(keep, function(keep) {
                return(data[keep,]) # mcols are carried over
            })

            ## merge nearby windows (within bin_size) to get broader TSS
            ## final fold change = avgFC of windows
            merged <- lapply(filtered.data, function(d) {
                dr <- GenomicRanges::granges(d)
                dr_reduced <- GenomicRanges::reduce(dr,
                                                    min.gapwidth = mergeLength,
                                                    ignore.strand = FALSE,
                                                    with.revmap = TRUE)

                mcols(dr_reduced) <- aggregate(dr,
                                               mcols(dr_reduced)$revmap,
                                               score = BiocGenerics::mean(filter))
                return(dr_reduced)
            })

            # update the Capset object
            names(merged) <- unique(as.character(groups))
            CSobject@tss_detected <- GenomicRanges::GRangesList(merged)

            ## Calculate prop reads in TSS per group
            message("Counting reads within detected TSS")
            mergedall <- base::Reduce(S4Vectors::union, merged)
            si$num_intss <- as.numeric(numReadsInBed(mergedall, bam.files, countall = countall))
            sampleInfo(CSobject) <- si

            # Add the results as a list and save as .Rdata
            output <- list(
                counts.windows = data,
                counts.background = wider)

            if (!(is.null(outfile_prefix))) {
                message("Writing filtering information as .Rdata")
                save(output, file = paste0(outfile_prefix, ".Rdata"))
            }

            return(CSobject)
    })

#' Export the detected TSS from CapSet object as .bed files
#'
#' @rdname exportTSS
#' @param CSobject The modified CapSet object after running \code{\link{detectTSS}} function
#' @param outfile_prefix Prefix (with path) for output .bed files
#' @param pergroup If TRUE, write output per group of samples
#' @param merged If TRUE, write merged bed file (union of all groups)
#'
#' @return .bed file(s) containing detected TSS.
#'
#' @importFrom rtracklayer export.bed
#' @export
#' @examples
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' # export tss
#' exportTSS(cs, merged = TRUE, outfile_prefix = "testTSS")
#'

setMethod("exportTSS",
          signature = "CapSet",
          function(CSobject,
                   outfile_prefix,
                   pergroup,
                   merged) {
              mergedBED <- CSobject@tss_detected
              if (isTRUE(pergroup)) {
                  ## write merged output for each group
                  message("Writing output .bed files per group")
                  mapply(
                      function(bedfile, group) {
                          export.bed(object = bedfile, con = group)
                      },
                      bedfile = mergedBED,
                      group = paste0(outfile_prefix, "_" , names(mergedBED), ".bed")
                  )

              }
              if (isTRUE(merged)) {
                  ## write out the union of GRanges
                  message("Writing merged .bed files")
                  mergedall <- base::Reduce(S4Vectors::union, mergedBED)
                  export.bed(mergedall,
                            con = paste(outfile_prefix, "merged.bed", sep = "_"))
              }

          })
