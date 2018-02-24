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
