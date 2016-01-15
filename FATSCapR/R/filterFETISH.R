#### ~~~~ Part of the FATSCapR package for analysis of FETISH data ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)


#' Split the composite BAM file using internal indexes
#'
#' @param bamFile Input BAM file
#' @param index_list A list of indexes for splitting (character vector)
#' @param outfile_list A list of output file names (with order corresponding to that of index_list)
#' @param rep_indexList A list of replicate indexes (if you want to split by the replicates too)
#' @param nthreads Number of threads
#'
#' @return Filtured files
#' @export
#'
#' @examples
#' filter_BAM(bamFile = "FATSCapR/tests/testthat/test_sorted.bam",
#'  	   index_list = c("TTAGCC" ,"CAAGTG"),
#'  	    	   outfile_list = c("FATSCapR/tests/testthat/test_filt1.bam",
#'  	    	   "FATSCapR/tests/testthat/test_filt2.bam"), nthreads = 10)


split_BAM <- function(bamFile, index_list, outfile_list, rep_indexList = NULL, nthreads) {

	## Write a CLOSURE that return the function to search idx in readname
	message("Creating Filtering Rules")
	make_FilterFunc <- function(idx_name,rep_name = NULL){

		function(df){
			df$qname <- as.character(df$qname)
			df_sep <- data.frame(idx = vapply(strsplit(df$qname, "#"), "[[", character(1), 2),
					     stringsAsFactors = FALSE)
			df_sep2 <- data.frame(idx = vapply(strsplit(df_sep$idx, ":"), "[[", character(1), 1),
					      stringsAsFactors = FALSE)
#			df_sep3 <- data.frame(idx = vapply(strsplit(df_sep$idx, ":"), "[[", character(1), 3),
#					      stringsAsFactors = FALSE)

			## grep for idx and repIdx and combine them (AND)
			# If the replicate index name is given, then grep that from readName and
			# return TRUE if both idx and repIdx is present. Else, only check for idx.

#			grep_idx_res <- grepl(idx_name,df_sep2$idx)
#			grep_repidx_res = ifelse( rep_name != NULL,
#						  grepl(rep_name,df_sep3$idx), TRUE )


#			return(grep_idx_res & grep_repidx_res)
			return(grepl(idx_name,df_sep2$idx))
		}

	}

	## A function to create filterRules argument
	make_FilterRules <- function(FilterFunc){
			return(FilterRules(list(FilterFunc)) )
	}


	## Now put them together, get lists back
	# I have yet not implemented the version where it works for both idx and rep_idx.
	# have to do it (probably using mapply)
	filtfuncs <- lapply(index_list,make_FilterFunc)
	filtrules <- lapply(filtfuncs, make_FilterRules)

	## Filter the files in parallel
	message("Filtering the BAM file")
	destinations <- outfile_list

	register(MulticoreParam(workers = nthreads))

	bplapply(seq_along(destinations), function(i, file, destinations, filtrules) {
		filterBam(file, destinations[i], filter = filtrules[[i]])
	}, bamFile, destinations, filtrules)
	## Files written
	message("Done!")
}
