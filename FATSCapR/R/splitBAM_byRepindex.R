#### ~~~~ Part of the FATSCapR package for analysis of FETISH data ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)


#' Split the composite BAM file using replicate indexes
#'
#' @param bamFile Input BAM file
#' @param repindex_list A list of replicate indexes for splitting (character vector)
#' @param outfile_list A list of output file names (with order corresponding to that of repindex_list)
#' @param nthreads Number of threads
#'
#' @return Filtered files by replicate Index
#' @export
#'
#' @examples
#' 
#' splitBAM_byRepindex(bamFile = "FATSCapR/tests/testthat/test_sorted.bam",
#'  	   repindex_list = c("AA" ,"AT"),
#'  	    	   outfile_list = c("FATSCapR/tests/testthat/test_filt3.bam",
#'  	    	   "FATSCapR/tests/testthat/test_filt4.bam"), nthreads = 10)



splitBAM_byRepindex <- function(bamFile, repindex_list, outfile_list, nthreads = 1) {
	
	## Write a CLOSURE that return the function to search idx in readname
	message("Creating Filtering Rules")
	make_FilterFunc <- function(rep_name){
		
		function(df){
			df$qname <- as.character(df$qname)
			df_sep <- data.frame(idx = vapply(strsplit(df$qname, "#"), "[[", character(1), 2),
						   stringsAsFactors = FALSE)
			df_sep3 <- data.frame(idx = vapply(strsplit(df_sep$idx, ":"), "[[", character(1), 3),
							stringsAsFactors = FALSE)

			return(grepl(rep_name,df_sep3))
		}
		
	}
	
	## A function to create filterRules argument
	make_FilterRules <- function(FilterFunc){
		return(FilterRules(list(FilterFunc)) )
	}
	
	
	## Now put them together, get lists back
	
	filtfuncs <- lapply(repindex_list,make_FilterFunc)
	filtrules <- lapply(filtfuncs, make_FilterRules)

	## Filter the files in parallel
	message("Filtering the BAM file")
	destinations <- ifelse(rep_indexList != NULL, paste(outfile_list, rep_indexList,sep = "_"), outfile_list)
	
	param = BiocParallel::MulticoreParam(workers = nthreads)
	bplapply(seq_along(destinations), function(i, file, destinations, filtrules) {
		filterBam(file, destinations[i], filter = filtrules[[i]])
	}, bamFile, destinations, filtrules,
	BPPARAM = param)
	
	## Files written
	message("Done!")
}
