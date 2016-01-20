#### ~~~~ Part of the FATSCapR package for analysis of FETISH data ~~~~ ####
### (c) Vivek Bhardwaj (bhardwaj@ie-freiburg.mpg.de)


#' Split the composite BAM file using replicate indexes
#'
#' @param bamFile Input BAM file
#' @param outfile_prefix A prefix for output file (replicates IDs will be added as RR/YY)
#' @param nthreads Number of threads
#'
#' @return Filtered files by replicate Index
#' @export
#'
#' @examples
#' 
#' splitBAM_byRepindex(bamFile = "FATSCapR/tests/testthat/test_sorted.bam",
#'  	   	    	   outfile_prefix = "TEST", nthreads = 10)
#'  	   	    	   



splitBAM_byRepindex <- function(bamFile, outfile_prefix, nthreads = 1) {
	
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
	
	## Splitting by PURINES (RR) or PYRIMIDINES (YY)
	repindex_list <- list(RR = paste("AA","GG","GA","AG", sep = "|"),
				    YY = paste("CC","TT","CT","TC", sep = "|"))
	## Now put them together, get lists back
	
	filtfuncs <- lapply(repindex_list,make_FilterFunc)
	filtrules <- lapply(filtfuncs, make_FilterRules)

	## Filter the files in parallel
	message("Filtering the BAM file")
	destinations <- paste(outfile_prefix, names(repindex_list),sep = "_")
	
	param = BiocParallel::MulticoreParam(workers = nthreads)
	BiocParallel::bplapply(seq_along(destinations), function(i, file, destinations, filtrules) {
		filterBam(file, destinations[i], filter = filtrules[[i]])
	}, bamFile, destinations, filtrules,
	BPPARAM = param)
	
	## Files written
	message("Done!")
}
