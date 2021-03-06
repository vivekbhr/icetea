
#' @rdname sampleInfo
#' @export
setGeneric("sampleInfo", function(object, ...)
    standardGeneric("sampleInfo"))

#' @rdname sampleInfo
#' @export
setGeneric("sampleInfo<-", function(object, ..., value)
    standardGeneric("sampleInfo<-"))

#' Retrieve and replace sample information of a CapSet object
#'
#' @name sampleInfo
#' @param object The \code{\link{CapSet}} object
#' @param value Replacement DataFrame object
#' @param ... Additional options
#'
#' @return sample information data.frame
#' @export
#' @examples
#'
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' # get sampleinfo
#' si <- sampleInfo(cs)
#' # modify
#' si$samples <- paste0("sample_", seq_along(1:nrow(si)) )
#' # replace
#' sampleInfo(cs) <- si
NULL
#> NULL

#' @rdname plotTSSprecision
#' @export
setGeneric("plotTSSprecision",
            function(reference,
                    detectedTSS,
                    distanceCutoff = 500,
                    outFile = NULL,
                    ...) {
                standardGeneric("plotTSSprecision")
            })

#' @name plotReadStats
#' @rdname plotReadStats
#' @export
setGeneric("plotReadStats",
           function(CSobject,
                    plotType = "dodge",
                    plotValue = "numbers",
                    outFile = NULL) {
                standardGeneric("plotReadStats")
            })

#' @name demultiplexFASTQ
#' @rdname demultiplexFASTQ
#' @export
setGeneric("demultiplexFASTQ",
           function(CSobject,
                    outdir,
                    max_mismatch = 0,
                    ncores = 1)
            standardGeneric("demultiplexFASTQ"))

#' @name mapCaps
#' @rdname mapCaps
#' @export
setGeneric("mapCaps",
           function(CSobject,
                    genomeIndex,
                    outdir,
                    externalGTF = NULL,
                    ncores = 1,
                    logfile = NULL)
            standardGeneric("mapCaps"))

#' @name filterDuplicates
#' @rdname filterDuplicates
#' @export
setGeneric("filterDuplicates",
           function(CSobject,
                    outdir,
                    ncores = 1,
                    keepPairs = FALSE)
    standardGeneric("filterDuplicates"))

#' @name detectTSS
#' @rdname detectTSS
#' @export
setGeneric("detectTSS",
           function(CSobject,
                    groups,
                    outfile_prefix = NULL,
                    windowSize = 10L,
                    sliding = TRUE,
                    foldChange = 2,
                    mergeLength = 1L,
                    restrictChr = NULL,
                    ncores = 1,
                    readPos = "start")
            standardGeneric("detectTSS"))

#' @name exportTSS
#' @rdname exportTSS
setGeneric("exportTSS",
           function(CSobject,
                    outfile_prefix,
                    pergroup = FALSE,
                    merged = TRUE)
            standardGeneric("exportTSS"))

#' @name getGeneCounts
#' @rdname getGeneCounts
#' @export
setGeneric("getGeneCounts",
            function(CSobject,
                transcriptGRL,
                regionAroundTSS = 500,
                outfile = NA,
                ncores = 1)
            standardGeneric("getGeneCounts"))

#' @name getNormFactors
#' @rdname getNormFactors
#' @export
setGeneric("getNormFactors",
            function(CSobject,
                    features,
                    method = "TMM",
                    ...)
            standardGeneric("getNormFactors"))

#' Detect differentially expressed Transcription Start Sites between two conditions (fit model)
#' @name fitDiffTSS
#' @rdname fitDiffTSS
#' @export
setGeneric("fitDiffTSS",
            function(CSobject,
                    TSSfile = NULL,
                    groups,
                    method = "DESeq2",
                    normalization = NULL,
                    normFactors = NULL,
                    outplots = NULL,
                    plotRefSample = NA,
                    ncores = 1)
            standardGeneric("fitDiffTSS"))

#' Detect differentially expressed Transcription Start Sites between two conditions (test)
#' @name detectDiffTSS
#' @rdname detectDiffTSS
#' @param fit DGEGLM object (output of \code{\link{fitDiffTSS}} command )
#' @param testGroup Test group name
#' @param contGroup Control group name
#' @param TSSfile The TSS .bed file used for \code{\link{fitDiffTSS}} command (if method "edgeR" was used)
#' @param MAplot_fdr FDR threshold to mark differentially expressed TSS in MAplot (NA = Don't make an MAplot)
#'
#' @return A \code{\link{GRanges}} object containing p-values of differential expression for each TSS.
#' @export
setGeneric("detectDiffTSS",
           function(fit,
                    testGroup,
                    contGroup,
                    TSSfile = NULL,
                    MAplot_fdr = NA)
            standardGeneric("detectDiffTSS"))
