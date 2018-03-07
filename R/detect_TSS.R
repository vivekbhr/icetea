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


#' Detection of Trancription start sites based on local enrichment
#'
#' @rdname detectTSS
#' @param CSobject CapSet object created using \code{\link{newCapSet}} function
#' @param groups a character vector that contains group name of the sample, for replicate-based TSS
#'               calling (see example)
#' @param outfile_prefix Output name prefix for the .Rdata file containing window counts, background counts
#'                       and filtering statistics calculated during TSS detection.
#' @param foldChange A fold change cutoff of local enrichment to detect the TSS. For samples with
#'        usual' amount of starting material and squencing depth (>=5ug starting material,
#'        = 5 mil reads/sample), a cut-off of 6 fold can be used. For samples with low
#'        amount of material or sequencing depth, use a lower cut-off (eg. use 2-fold for
#'        samples with 500ng starting material).
#' @param restrictChr Chromosomes to restrict the analysis to.
#'
#'
#' @return .bed files containing TSS position for each group, along with a bed file for consensus
#'        (union) TSS sites of all samples.
#'
#' @export
#' @importFrom utils write.table
#' @importFrom SummarizedExperiment colData rowRanges
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
#'            foldChange = 6, restrictChr = "X")
#'

setMethod("detectTSS",
          signature = "CapSet",
          function(CSobject,
                   groups,
                   outfile_prefix,
                   foldChange,
                   restrictChr) {
              # check whether group and outfile_prefix is provided
              if (missing(outfile_prefix))
                  stop("Please provide outfile_prefix!")
              if (missing(groups))
                  stop("Please provide groups!")

              # convert group to char
              si <- sampleInfo(CSobject)
              design <-
                  data.frame(row.names = si$samples, group = as.character(groups))

              if (is.null(si$filtered_file)) {
                  message("Filtered files not found under sampleInfo(CSobject). Using mapped files")
                  bam.files <- si$mapped_file
              } else {
                  bam.files <- si$filtered_file
              }
              if (sum(file.exists(bam.files)) != length(bam.files)) {
                  stop("One or more bam files don't exist! Check sampleInfo(CSobject) ")
              }

              # Define read params
              frag.len <- NA
              win.width <- 10
              param <-
                  readParam(minq = 2,
                            forward = NULL,
                            restrict = restrictChr)
              regionparam <- readParam(minq = 2, restrict = restrictChr)

              # Count reads into sliding windows
              data <-
                  strandedCounts(
                      bam.files,
                      param = param,
                      ext = frag.len,
                      width = win.width,
                      bin = TRUE
                  )
              colnames(data) <- rownames(design)
              colData(data) <- c(colData(data), design)

              # Get counts for 2kb local region surrounding each bin
              surrounds <- 2000
              neighbor <- suppressWarnings(GenomicRanges::trim(
                  GenomicRanges::resize(rowRanges(data),
                                        surrounds, fix = "center")
              ))

              wider <- suppressWarnings(regionCounts(
                  bam.files,
                  param = regionparam,
                  regions = neighbor,
                  ext = frag.len
              ))

              colnames(wider) <- rownames(design)
              colData(wider) <- c(colData(wider), design)

              ## take out groups --> Generate filter statistics for each group (based on local enrichment)
              filterstat <- lapply(unique(design$group), function(x) {
                  stat <- filterWindows(data[, data$group == x],
                                        wider[, wider$group == x],
                                        type = "local")
                  return(stat)
              })

              # Require X-fold enrichment over local background to keep the window (similar to MACS)
              keep <- lapply(filterstat, function(x) {
                  kp <- x$filter > log2(foldChange)
                  return(kp)
              })

              filtered.data <- lapply(keep, function(keep) {
                  return(data[keep,])
              })

              ## merge nearby windows (within 10bp) to get broader TSS
              merged <- lapply(filtered.data, function(d) {
                  return(mergeWindows(d, tol = 10L, ignore.strand = FALSE))
              })
              # update the Capset object
              merged <- lapply(merged, function(x)
                  return(x$region))
              names(merged) <- unique(as.character(groups))
              CSobject@tss_detected <- GenomicRanges::GRangesList(merged)

              ## Calculate prop reads in TSS per group
              message("Counting reads within detected TSS")
              mergedall <- base::Reduce(S4Vectors::union, merged)
              si$num_intss <- as.numeric(numReadsInBed(mergedall, bam.files))
              sampleInfo(CSobject) <- si

              # Add the results as a list and save as .Rdata
              output <- list(
                  counts.windows = data,
                  counts.background = wider,
                  filter.stats = S4Vectors::DataFrame(filterstat[[1]])
              )
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
                          rtracklayer::export.bed(object = bedfile, con = group)
                      },
                      bedfile = mergedBED,
                      group = paste0(outfile_prefix, "_" , names(mergedBED), ".bed")
                  )

              }
              if (isTRUE(merged)) {
                  ## write out the union of GRanges
                  message("Writing merged .bed files")
                  mergedall <- base::Reduce(S4Vectors::union, mergedBED)
                  rtracklayer::export.bed(mergedall,
                                          con = paste(outfile_prefix, "merged.bed", sep = "_"))
              }

          })
