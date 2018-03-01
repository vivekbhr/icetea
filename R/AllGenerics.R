#' @rdname plot_TSSprecision
#' @export
setGeneric("plot_TSSprecision",
       function(reference, detectedTSS,
           distanceCutoff = 500, outFile = NULL, ...) {
           standardGeneric("plot_TSSprecision")
       })

#' @rdname sampleInfo
#' @export
setGeneric("sampleInfo", function(object,...) standardGeneric("sampleInfo"))

#' @rdname sampleInfo
#' @export
setGeneric("sampleInfo<-", function(object,...,value) standardGeneric("sampleInfo<-"))

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
