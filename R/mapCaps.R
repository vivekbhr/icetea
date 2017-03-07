#### ~~~~ Part of the mapcapR package for analysis of FETISH data ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)

#' Map the paired-end MAPCap data.
#'
#' @param index character string giving the basename of the Subread index file.
#' @param R1 forward read (_R1) fastq, output of \code{\link{trimFastqIndex}} command.
#' @param R2 reverse read (_R2) fastq, output of \code{\link{trimFastqIndex}} command.
#' @param output output file name (without ".bam" extention)
#' @param nthreads number of threads to use for mapping.
#' @param logfile a log file to write the processing message.
#' @param ... additional arguments passed to the RSubread::align function.
#' @return bamfile A mapped BAM file for the sample.
#'
#' @examples
#' r1 <- system.file("extdata", "testout_R1.fastq.gz", package = "mapcapR")
#' r2 <- system.file("extdata", "testout_R2.fastq.gz", package = "mapcapR")
#' \dontrun{
#' mapCaps(index,R1 = r1, R2 = r2, output = "test_mapped", nthreads = 10, logfile=NULL)
#' }
#'
#' @export
#'

mapCaps <- function(index, R1, R2, output, nthreads, logfile = NULL,...){
	# open a logfile if given
	if(!is.null(logfile)){
		sink(logfile)
	}

	# test for trimmed R2
	cat("Checking for trimmed R2\n\n")
	read2 <- gzfile(R2)
	data <- readLines(read2,100)
	close(read2)

	# Check if the read header has "#" (which is introduced during trimming.)
	header <- data[seq(1,100,4)]
	if (unique(grepl("#",header)) != TRUE) {
		stop("Stop! read R2 seems untrimmed. Run trimFETISH first.")
	}

	cat("Mapping the FETISH data\n\n")
	# Align using RSubread
	tmpout <- paste0(output,".tmp.bam")
	Rsubread::align(index = index,
			readfile1 = R1,
			readfile2 = R2,
			output_file = tmpout,
			nthreads = nthreads,
			...)

	# Sort and Index
	cat("Sorting and Indexing")
	Rsamtools::sortBam(file = tmpout, destination = output)# adds .bam suffix
	Rsamtools::index(paste0(output, ".bam") )
	file.remove(tmpout)

	# Close logfile
	if(!is.null(logfile)){
		sink()
	}
}
