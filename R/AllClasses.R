#' Create a new CapSet object
#'
#' @param expMethod experiment method ('CAGE', 'RAMPAGE' or 'MAPCap')
#' @param fastqType fastq file type. Should be either 'single' or 'paired'
#' @param fastq_R1 path for Read R1 (or file path for single end reads)
#' @param fastq_R2 path for Read R2 (for paired end reads)
#' @param sampleInfo data.frame with sample information, with row.names = barcode sequences,
#'                   and first column containing corresponding sample names.
#'
#' @return An object of class CapSet
#' @importFrom methods new
#' @export
#'
#' @examples
#'
#' # list of barcode IDs
#' idxlist <- c("CAAGTG", "AGATGC", "TGTGAG",
#' 		 "GGTTAC", "TTAGCC", "AGTCGA")
#' # corresponding sample names
#' fnames <-  c("WTa", "WTb", "WTc",
#' 		 "MLEa", "MLEb", "MLEc")
#'
#' # create a new capset object and save
#' dir <- system.file("extdata", package="icetea")
#'
#' cs <- newCapSet(expMethod = 'MAPCap', fastqType = 'paired',
#'		fastq_R1 = file.path(dir, 'mapcap_test_R1.fastq.gz'),
#'		fastq_R2 = file.path(dir, 'mapcap_test_R2.fastq.gz'),
#'		sampleInfo = data.frame(row.names = idxlist, samples = fnames))
#'
#' save(cs, file = file.path(dir, "CSobject.Rdata") )
#'

newCapSet <- function(sampleInfo, expMethod, fastqType, fastq_R1, fastq_R2 = NULL) {
    stopifnot(class(sampleInfo) %in% c("data.frame", "DataFrame"))
    # rename columns
    colnames(sampleInfo) <- "samples"
    # convert sampleInfo to a DataFrame
    info <- S4Vectors::DataFrame(sampleInfo)
    # create an instance of CapSet
    new("CapSet",
        sampleInfo = info,  # sample Information
        fastqType = fastqType,
        fastq_R1 = fastq_R1,
        fastq_R2 = fastq_R2,
        expMethod = expMethod,
        tss_detected = NULL)
}

#' Check capset validity
#'
#' @param object capset object
#'
#' @return boolean
#'

check_capSet <- function(object) {
    errors <- character()

    ## extract slots
    fqtype <- object@fastqType
    R1 <- object@fastq_R1
    R2 <- object@fastq_R2
    exp <- object@expMethod
    info <- object@sampleInfo
    tss <- object@tss_detected

    ## validate slots
    # fastq
    if(!(fqtype %in% c("single", "paired") )) {
    msg <- paste0("Wrong fastq type : ", fqtype, ". Should be either 'single' or 'paired' ")
    errors <- c(errors, msg)

    } else if(fqtype == "single" & !file.exists(R1) ) {
    msg <- paste0("Please specify correct fastq file path for fastq_R1 ")
    errors <- c(errors, msg)

    } else if(fqtype == "paired" & !(file.exists(R1) | file.exists(R2) ) ) {
    msg <- paste0("Please specify correct fastq file path for both fastq_R1 and fastq_R2 ")
    errors <- c(errors, msg)
    }
    # experiment
    if(!(exp %in% c("CAGE", "RAMPAGE", "MAPCap") )) {
    msg <- paste0("Experiment type should be among : 'CAGE', 'RAMPAGE' or 'MAPCap' ")
    errors <- c(errors, msg)

    }
    # sampleInfo
    if(!(class(info) %in% c("data.frame", "DataFrame")) ) {
    msg <- paste0("sampleInfo should be a data frame ")
    errors <- c(errors, msg)
    }
        # TSS info
        if(!(class(tss) %in% c("NULL", "GRangesList")) ) {
        	msg <- paste0("tss_detected should be a GRangesList object ")
        	errors <- c(errors, msg)
        }

    ## return
    if (length(errors) == 0) TRUE else errors
}

## char or NULL class
setClassUnion("charOrNULL", c("character", "NULL"))

#' CapSet object
#'
#' @rdname newCapSet
#' @importClassesFrom S4Vectors DataFrame
#'
CapSet <- setClass("CapSet",
       slots = c(fastqType = "character",
        fastq_R1 = "character",
        fastq_R2 = "charOrNULL",
        expMethod = "character",
        sampleInfo = "DataFrame",
        tss_detected = "ANY"
        ),
       validity = check_capSet)
