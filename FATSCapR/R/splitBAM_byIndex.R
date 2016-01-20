#### ~~~~ Part of the FATSCapR package for analysis of FETISH data ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)


#' Split the composite BAM file using internal indexes
#'
#' @param bamFile Input BAM file
#' @param index_list A list of indexes for splitting (character vector)
#' @param outfile_list A list of output file names (with order corresponding to that of index_list)
#' @param max_mismatch No. of mismatches allowed in index (maxium 1 recommended)
#' @param nthreads Number of threads
#'
#' @return Filtered files
#' @export
#'
#' @examples
#' splitBAM_byIndex(bamFile = "FATSCapR/tests/testthat/test_sorted.bam",
#'  	   index_list = c("TTAGCC" ,"CAAGTG"),
#'  	    	   outfile_list = c("FATSCapR/tests/testthat/test_filt1.bam",
#'  	    	   "FATSCapR/tests/testthat/test_filt2.bam"), nthreads = 10)


splitBAM_byIndex <- function(bamFile, index_list, outfile_list, max_mismatch = 0, nthreads = 1) {

	## Write a CLOSURE that return the function to search idx in readname
	message("Creating Filtering Rules")
	make_FilterFunc <- function(idx_name,maxM = max_mismatch){
		
		idx_name = Biostrings::DNAString(idx_name)
		
		function(df){
			df$qname <- as.character(df$qname)
			df_sep <- data.frame(idx = vapply(strsplit(df$qname, "#"), "[[", character(1), 2),
					     stringsAsFactors = FALSE)
			df_sep2 <- data.frame(idx = vapply(strsplit(df_sep$idx, ":"), "[[", character(1), 1),
					      stringsAsFactors = FALSE)
			df_sep2 <- Biostrings::DNAStringSet(df_sep2$idx)

			grep_idx_res <- as.logical(Biostrings::vcountPattern(idx_name,df_sep2,max.mismatch = maxM))

			return(grep_idx_res)
		}

	}

	## A function to create filterRules argument
	make_FilterRules <- function(FilterFunc){
			return(S4Vectors::FilterRules(list(FilterFunc)))
	}


	## Now put them together, get lists back

	filtfuncs <- lapply(index_list,make_FilterFunc)
	filtrules <- lapply(filtfuncs, make_FilterRules)

	## Filter the files in parallel
	message("Filtering the BAM file")
	destinations <- outfile_list
	
	param = BiocParallel::MulticoreParam(workers = nthreads)
	BiocParallel::bplapply(seq_along(destinations), function(i, file, destinations, filtrules) {
								filterBam(file, destinations[i], filter = filtrules[[i]])
								}, bamFile, destinations, filtrules,
		   BPPARAM = param)
	
	## Files written
	message("Done!")
}
