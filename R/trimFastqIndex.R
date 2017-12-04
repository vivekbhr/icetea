#' Trim the sample indexes off the fastq files (RAMPAGE/MAPCap)
#'
#' @param CapSet The CapSet object
#' @param outdir output directory to keep the trimmed fastq files
#'
#' @return fastq files
#' @export
#'
#' @examples
#'
trimFastqIndex <- function(CapSet, outdir){

	# check if output dir exists
	if(!dir.exists(outdir)) dir.create(outdir)
	outfile_name <- "trimmed"

	input_R1 <- CapSet@fastq_R1
	input_R2 <- CapSet@fastq_R2
	trimmed_R1 <- file.path(outdir, "trimmed_R1.fastq.gz")
	trimmed_R2 <- file.path(outdir, "trimmed_R2.fastq.gz")

	message("Trimming the barcodes : ")
	# passes the arguments to the C function, which accepts two chars as double pointers
	if (CapSet@expMethod == "MAPCap") {
		message("MAPCap")
		.C("trimFastq", input_R1, input_R2, trimmed_R1, trimmed_R2)
	} else if (CapSet@expMethod == "RAMPAGE") {
		message("RAMPAGE")
		.C("trimFastq_rampage", input_R1, input_R2, trimmed_R1, trimmed_R2)
	} else if (CapSet@expMethod == "CAGE") {
		warning("No barcode positions known for the conventional CAGE protocol. No trimming performed ")
		trimmed_R1 <- input_R1
		trimmed_R2 <- input_R1
	}

	# output the modified CapSet object
	CapSet@trimmed_R1 <- trimmed_R1
	CapSet@trimmed_R2 <- trimmed_R2
	validObject(CapSet)
	return(CapSet)
}
