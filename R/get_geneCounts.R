

#' Get gene-level counts from TSS data
#'
#' @param txdb A txdb object
#' @param bamfiles Bam files to count the reads from
#' @param single_end whether reads are single end
#' @param outfile Tab-separated output file name (if required)
#'
#' @return data.frame with gene-level counts for all genes in the txdb object
#' @export
#'
#' @examples
#'

get_geneCounts <- function(txdb, bamfiles, single_end = TRUE,outfile = NA) {


	# get 500bp region around transcripts to count the reads
	dm6trans <- GenomicFeatures::transcriptsBy(txdb, "gene")

	#add gene names and unlist
	dm6trans <- mapply(function(x, name) {
		x$geneID <- name
		return(x)
		}, dm6trans, names(dm6trans))

	dm6trans <- unlist(GenomicRanges::GRangesList(dm6trans))
	dm6trans_resized <- GenomicRanges::resize(dm6trans, 500, fix = "start")

	# count reads
	tsscounts <- GenomicAlignments::summarizeOverlaps(dm6trans_resized,
							  reads = Rsamtools::BamFileList(bamfiles),
							  singleEnd = single_end)
	tsscounts.df <- SummarizedExperiment::assay(tsscounts)

	# sum all the TSS counts to gene counts
	id <- rownames(tsscounts.df)
	tsscounts.gene <- lapply(split(as.data.frame(tsscounts.df), id), colSums)
	tsscounts.gene <- do.call(rbind, tsscounts.gene)
	rm(tsscounts.df)

	# write back
	if(!is.na(outfile)) {
		write.table(tsscounts.gene,
			    file = outfile,
			    sep = "\t", quote = FALSE)
	}
	return(tsscounts.gene)

}
