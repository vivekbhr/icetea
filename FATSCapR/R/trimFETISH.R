#### ~~~~ Part of the FATSCapR package for analysis of FETISH data ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)

#' Trim the R2 of FETISH data and tag the read headers
#'
#' @param infile Input _R2.fastq (or _R2.fastq.gz) file.
#' @param outfile Name of the Output fastq file.
#' @return trimmedFile : The trimmed R2.fastq file.
#' @examples
#' trimFETISH("testdata/test_R2.fastq.gz","test_trimmed_R2.fastq")
#'

trimFETISH <- function(infile,outfile){
	message("Trimming the FETISH data")
	# passes the arguments to the C function, which accepts two chars as double pointers
	out <- .C("trimFETISH",as.character(infile),as.character(outfile))
	# output is a list with input filenames
	names(out) <- c("Input","Output")
	return(out) # Just print names on the screen
}
