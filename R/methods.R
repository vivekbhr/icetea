## show method
setMethod("show", "CapSet", function(object) {
    cat("An object of class CapSet", "\n")
    cat("---------------------------", "\n\n")
    cat("Experiment method : ", object@expMethod, "\n")
    #cat("FASTQ Type : ", object@fastqType, "\n")
    cat("FASTQ Read 1 : ", object@fastq_R1, "\n")
    cat("FASTQ Read 2 : ", object@fastq_R2, "\n")
    cat("\n", "Sample information : ", "\n")
    cat("-------------------------", "\n")
    print(object@sampleInfo)

    cat("\n", "TSS enrichment information : ", "\n")
    cat("-----------------------------", "\n")
    cat("Detected TSS per group", "\n")
    print(object@tss_detected)
})



#' @rdname sampleInfo
#' @export
#'
setMethod("sampleInfo",
          signature = "CapSet",
          function(object) {
              return(object@sampleInfo)
          })

#' @rdname sampleInfo
#' @exportMethod "sampleInfo<-"
#'
setReplaceMethod("sampleInfo",
                 signature = "CapSet",
                 function(object, value) {

                    df <- S4Vectors::DataFrame(value)
                    object@sampleInfo <- df
                    validObject(object)
                    return(object)
                 })
