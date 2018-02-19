#' @rdname plot_TSSprecision
#' @export
setGeneric("plot_TSSprecision",
	   function(reference, detectedTSS,
	   	 distanceCutoff = 500, outFile = NULL, ...) {
	   	standardGeneric("plot_TSSprecision")
	   })
