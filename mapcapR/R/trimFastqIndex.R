#### ~~~~ Part of the mapcapR package for analysis of FETISH data ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)

#' Trim the raw fastq files and tag the read headers
#'
#' @param input_R1 Input _R1.fastq (or _R1.fastq.gz) file.
#' @param input_R2 Input _R2.fastq (or _R2.fastq.gz) file.
#' @param output_R1 Name of the Output _R1.fastq.gz file.
#' @param output_R2 Name of the Output _R2.fastq.gz file.
#' @return output_R1, output_R2 : The trimmed R1 and R2 .fastq.gz files.
#' @examples
#' r1 <- system.file("extdata", "test_R1.fastq.gz", package = "mapcapR")
#' r2 <- system.file("extdata", "test_R2.fastq.gz", package = "mapcapR")
#' trimFastqIndex(r1, r2,"test_trimmed_R1.fastq.gz","test_trimmed_R2.fastq.gz")
#'
#' @useDynLib mapcapR trimFastq
#' @export

trimFastqIndex <- function(input_R1,input_R2,output_R1,output_R2){
	message("Trimming the mapcap data")
	# passes the arguments to the C function, which accepts two chars as double pointers
	out <- .C("trimFastq",as.character(input_R1),as.character(input_R2),
		  as.character(output_R1),as.character(output_R2))
	# output is a list with input filenames
	names(out) <- c("Input_R1","Input_R2","Output_R1","Output_R2")
	return(out) # Just print names on the screen
}
