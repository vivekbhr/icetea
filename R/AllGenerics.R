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
             nthreads = 1,
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
                    foldChange = 2,
                    restrictChr = NULL,
                    ncores = 1)
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
                single_end = TRUE,
                outfile = NA)
            standardGeneric("getGeneCounts"))

#' @name fitDiffTSS
#' @rdname fitDiffTSS
#' @export
setGeneric("fitDiffTSS",
            function(CSobject,
                    TSSfile = NULL,
                    groups,
                    normalization = "internal",
                    CSobjectSpikeIn = NULL,
                    outplots = NULL,
                    plotref)
            standardGeneric("fitDiffTSS"))

#' @name detectDiffTSS
#' @rdname detectDiffTSS
#' @export
setGeneric("detectDiffTSS",
           function(fit,
            testGroup,
            contGroup,
            TSSfile,
            MAplot_fdr = NA)
            standardGeneric("detectDiffTSS"))
