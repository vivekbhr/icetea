#### ~~~~ Part of the mapcapR package for analysis of MAPCap data ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)

#' Map the paired-end MAPCap data.
#'
#' @param index character string giving the basename of the Subread index file.
#' @param R1 forward read (R1) fastq file, output of \code{\link{trimFastqIndex}} command.
#' @param R2 reverse read (R2) fastq file, output of \code{\link{trimFastqIndex}} command.
#' @param outprefix output file prefix (without ".bam" extention)
#' @param nthreads number of threads to use for mapping.
#' @param logfile a log file to write the processing message.
#' @param ... additional arguments passed to the \code{\link{subjunc}} function.
#' @return bamfile A mapped BAM file for the sample.
#'
#' @examples
#' r1 <- system.file("extdata", "testout_R1.fastq.gz", package = "mapcapR")
#' r2 <- system.file("extdata", "testout_R2.fastq.gz", package = "mapcapR")
#' \dontrun{
#' mapCaps(index,R1 = r1, R2 = r2, outprefix = "test_mapped", nthreads = 10, logfile=NULL)
#' }
#'
#' @export
#'

mapCaps <- function(index, R1, R2, outprefix, nthreads, logfile = NULL,...){
	# open a logfile if given
	if(!is.null(logfile)){
		sink(logfile)
	}

	# test for trimmed R2 index
	message("Checking for trimmed barcodes\n")
	read2 <- gzfile(R2)
	data <- readLines(read2,100)
	close(read2)

	# Check if the read header has "#" (which is introduced during trimming.)
	header <- data[seq(1,100,4)]
	if (unique(grepl("#",header)) != TRUE) {
		stop("Stop! read R2 seems untrimmed. Run trimFastqIndex first.")
	}

	message("Mapping the data\n")
	# Align using RSubread
	tmpout <- paste0(outprefix,".tmp.bam")
	Rsubread::subjunc(index = index,
			readfile1 = R1,
			readfile2 = R2,
			output_file = tmpout,
			nthreads = nthreads,
			minFragLength=10,
			reportAllJunctions = TRUE,
			...)

	# Sort and Index
	message("Sorting and Indexing")
	Rsamtools::sortBam(file = tmpout, destination = outprefix)# adds .bam suffix
	Rsamtools::indexBam(paste0(outprefix, ".bam") )
	file.remove(tmpout)

	# Close logfile
	if(!is.null(logfile)){
		sink()
	}
}
