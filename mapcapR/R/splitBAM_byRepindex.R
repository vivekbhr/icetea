#### ~~~~ Part of the mapcapR package for analysis of FETISH data ~~~~ ####
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
#' bam <- system.file("extdata", "test_mapped.bam", package = "mapcapR")
#' splitBAM_byRepindex(bamFile = bam, outfile_prefix = "TEST", nthreads = 10)
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

			return(grepl(rep_name,df_sep3$idx))
		}

	}

	## A function to create filterRules argument
	make_FilterRules <- function(FilterFunc){
		return(S4Vectors::FilterRules(list(FilterFunc)) )
	}

	## Splitting by PURINES (RR) or PYRIMIDINES (YY)
	repindex_list <- list(RR = paste("AA","GG","GA","AG", sep = "|"),
				    YY = paste("CC","TT","CT","TC", sep = "|"))
	## Now put them together, get lists back

	filtfuncs <- lapply(repindex_list,make_FilterFunc)
	filtrules <- lapply(filtfuncs, make_FilterRules)

	## Filter the files in parallel
	message("Filtering the BAM file")
	destinations <- paste0(outfile_prefix,"_", names(repindex_list),".bam")

	param = BiocParallel::MulticoreParam(workers = nthreads)
	BiocParallel::bplapply(seq_along(destinations), function(i, file, destinations, filtrules) {
							Rsamtools::filterBam(file, destinations[i], filter = filtrules[[i]])
							}, bamFile, destinations, filtrules,
				BPPARAM = param)

	## Files written
	message("Done!")
}
