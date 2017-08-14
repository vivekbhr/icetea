

#' Create a new CapSet object
#'
#' @param expMethod experiment method ('CAGE', 'RAMPAGE' or 'MAPCap')
#' @param fastqType fastq file type. Should be either 'single' or 'paired'
#' @param fastq_R1 path for Read R1 (or file path for single end reads)
#' @param fastq_R2 path for Read R2 (for paired end reads)
#' @param sampleBarcodes for de-multiplexing, a data.frame with two columns. First column should have the barcodes,
#'                       while the second column should have the corresponding sample names.
#'
#' @return An object of class CapSet
#' @export
#'
#' @examples
#'
#' \dontrun{
#' newCapSet(fastqType = "single", fastq_R1 = "out_R1.fastq.gz",
#' fastq_R2 = "out_R2.fastq.gz" , expMethod = "RAMPAGE", sampleBarcodes = df)
#'
#'}
#'

newCapSet <- function(expMethod, fastqType, fastq_R1, fastq_R2 = NA, sampleBarcodes) {
	# create an instance of CapSet
	new("capSet",
	    fastqType = fastqType,
	    fastq_R1 = fastq_R1,
	    fastq_R2 = fastq_R2,
	    expMethod = expMethod,
	    sampleBarcodes = sampleBarcodes)
}

#' Check capset validity
#'
#' @param object capset object
#'
#' @return boolean
#'
#' @examples
#'

check_capSet <- function(object) {
	errors <- character()

	## extract slots
	fqtype <- object@fastqType
	R1 <- object@fastq_R1
	R2 <- object@fastq_R2
	exp <- object@expMethod
	barcodeInfo <- object@sampleBarcodes

	## validate slots
	# fastq
	if(!(fqtype %in% c("single", "paired") )) {
		msg <- paste0("Wrong fastq type : ", fqtype, " . Should be either 'single' or 'paired' ")
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

	## return
	if (length(errors) == 0) TRUE else errors
}

CapSet <- setClass("CapSet",
		   slots = c(fastqType = "character",
		   	  fastq_R1 = "character",
		   	  fastq_R2 = "character",
		   	  expMethod = "character",
		   	  sampleBarcodes = "data.frame"),
		   validity = check_capSet)
