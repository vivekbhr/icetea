#### ~~~~ Part of the FATSCapR package for analysis of FETISH data ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)

#' Map the paired-end FETISH data.
#'
#' @param index character string giving the basename of the Subread index file.
#' Index files should be located in the current directory.
#' @param R1 forward read (_R1) fastq.
#' @param R2 reverse read (_R2) fastq.
#' @param output output file name
#' @param nthreads number of threads to use for mapping.
#' @param logfile a log file to write the processing message.
#' @param ... additional arguments passed to the RSubread::align function.
#' @return bamfile A mapped BAM file for the sample.
#' @examples
#'
#'

mapFETISH <- function(index,R1,R2,output,nthreads,logfile=NULL,...){
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
	Rsubread::align(index = index,
			readfile1 = R1,
			readfile2 = R2,
			output_file = output,
			nthreads = nthreads,
			...)
	# Close logfile
	if(!is.null(logfile)){
		sink()
	}
}
