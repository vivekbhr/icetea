## show method
setMethod("show", "CapSet", function(object) {
    cat("An object of class CapSet", "\n")
    cat("---------------------------", "\n\n")
    cat("Experiment method : ", object@expMethod, "\n")
    cat("FASTQ Type : ", object@fastqType, "\n")
    cat("FASTQ Read 1 : ", object@fastq_R1, "\n")
    cat("FASTQ Read 2 : ", object@fastq_R2, "\n")
    cat("\n", "Sample information : ", "\n")
    cat("-------------------------","\n")
    print(object@sampleInfo)

    cat("\n", "TSS enrichment information : ", "\n")
    cat("-----------------------------","\n")
    cat("Detected TSS per group", "\n")
    print(object@tss_detected)
} )

#' get sample information data frame
#' @param object the \code{\link{CapSet}} object
#' @param ... ...
#' @return sample information data.frame
#' @export
#' @examples
#'
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' # get sampleinfo
#' si <- sampleInfo(cs)
#'
setGeneric("sampleInfo", function(object,...) standardGeneric("sampleInfo"))

#' get sample information data frame
#'
#' @param object the \code{\link{CapSet}} object
#'
#' @docType methods
#' @return sample information data.frame
#' @export
#' @examples
#'
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' # get sampleinfo
#' si <- sampleInfo(cs)
#'

setMethod("sampleInfo",
      signature = "CapSet",
      function(object) {
      	return(object@sampleInfo)
      })

#' reset sample information data frame
#' @param object The \code{\link{CapSet}} object
#' @param ... ...
#' @param value new value
#' @export
#'
#' @return None
#' @examples
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' # get sampleinfo
#' si <- sampleInfo(cs)
#' # modify
#' si$samples <- paste0("sample_", 1:6)
#' # replace
#' sampleInfo(cs) <- si
#'
setGeneric("sampleInfo<-", function(object,...,value) standardGeneric("sampleInfo<-"))

#' reset sample information data frame
#' @param object The \code{\link{CapSet}} object
#' @param value The replacement sampleInfo data frame
#' @return None
#' @exportMethod "sampleInfo<-"
#'
#' @examples
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#' # get sampleinfo
#' si <- sampleInfo(cs)
#' # modify
#' si$samples <- paste0("sample_", 1:6)
#' # replace
#' sampleInfo(cs) <- si
#'
#'
setReplaceMethod("sampleInfo",
     signature = "CapSet", #value = c("data.frame", "DataFrame") ),
     function(object, value) {
     	df <- S4Vectors::DataFrame(value)
     	object@sampleInfo <- df
     	validObject(object)
     	return(object)
     })
