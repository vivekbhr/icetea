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
#' @export
#'
#' @examples
#'
#' \dontrun{
#' df <- data.frame(row.names = c("TTAGCC" ,"CAAGTG"), samples = c("one", "two") )
#' cs <- newCapSet(fastqType = "paired",
#'                 fastq_R1 = "inst/extdata/test_R1.fastq.gz",
#'                 fastq_R2 = "inst/extdata/test_R2.fastq.gz",
#'                 expMethod = "MAPCap", sampleInfo = df)
#'}
#'

newCapSet <- function(sampleInfo, expMethod, fastqType, fastq_R1, fastq_R2 = NULL) {
	stopifnot(!(class(sampleInfo) %in% c("data.frame", "DataFrame")))
	# convert sampleInfo to a DataFrame
	info <- S4Vectors::DataFrame(sampleInfo)
	# create an instance of CapSet
	new("CapSet",
	    sampleInfo = info,  # sample Information
	    fastqType = fastqType,
	    fastq_R1 = fastq_R1,
	    fastq_R2 = fastq_R2,
	    trimmed_R1 = NULL,
	    trimmed_R2 = NULL,
	    expMethod = expMethod,
	    counts.windows = NULL, # TSS : window counts
	    counts.background = NULL, # TSS : background counts
	    filter.stats = NULL) # TSS : filter stats

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

	wincounts <- object@counts.windows
	bgcounts <- object@counts.background
	filtstats <- object@filter.stats

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
	if (!(class(info) %in% c("data.frame", "DataFrame")) ) {
		msg <- paste0("sampleInfo should be a data frame ")
		errors <- c(errors, msg)
	}
        # TSS info
	if(!(class(wincounts) %in% c("RangedSummarizedExperiment", "NULL")) |
	   !(class(bgcounts) %in% c("RangedSummarizedExperiment", "NULL")) ) {
		msg <- paste0("window counts and background counts should either be or a
			      RangedSummarizedExperiment object")
		errors <- c(errors, msg)
	}

	if (!(class(filtstats) %in% c("NULL", "DataFrame")) ) {
		msg <- paste0("sampleInfo should be a data frame ")
		errors <- c(errors, msg)
	}

	## return
	if (length(errors) == 0) TRUE else errors
}

## char or NULL class
setClassUnion("charOrNULL", c("character", "NULL"))
## RSE or NULL class
#setClassUnion("RSEorNULL", c("RangedSummarizedExperiment", "NULL"))
## DataFrame or NULL class
#setClassUnion("DForNULL", c("DataFrame", "NULL"))

#' CapSet object
#'
#' @rdname newCapSet
#' @importClassesFrom S4Vectors DataFrame
#' @import SummarizedExperiment
#'
CapSet <- setClass("CapSet",
		   #contains = "DataFrame",
		   slots = c(fastqType = "character",
		   	  fastq_R1 = "character",
		   	  fastq_R2 = "charOrNULL",
		   	  trimmed_R1 = "charOrNULL",
		   	  trimmed_R2 = "charOrNULL",
		   	  expMethod = "character",
		   	  sampleInfo = "DataFrame",
		   	  counts.windows = "ANY", # TSS : window counts
		   	  counts.background = "ANY", # TSS : background counts
		   	  filter.stats = "ANY" # TSS : filter stats
		   	  ),
		   validity = check_capSet)
