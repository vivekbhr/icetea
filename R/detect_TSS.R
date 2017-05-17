

#' Detection of Trancription start sites based on local enrichment
#'
#' @param bam.files List of bam files to use
#' @param design A data frame with rownames = sample names and a column called 'group'
#' 		that contains information about the sample group (see example)
#'
#' @param foldChange A fold change cutoff of local enrichment to detect the TSS. For samples with
#' 		'usual' amount of starting material and squencing depth (>=5ug starting material,
#' 		>= 5 mil reads/sample), a cut-off of 6 fold can be used. For samples with low
#' 		amount of material or sequencing depth, use a lower cut-off (eg. use 2-fold for
#' 		samples with 500ng starting material).
#'
#' @param restrictChr Chromosomes to restrict the analysis to.
#'
#' @param outfile_prefix Output name prefix for the .bed files
#'
#'
#' @return .bed files containing TSS position for each group, along with a bed file for consensus
#' 	   (union) TSS sites of all samples.
#' @export
#'
#' @examples
#'
#'


detect_TSS <- function(bam.files, design, foldChange = 2, restrictChr, outfile_prefix) {

	# Define read params
	frag.len <- NA
	win.width <- 10
	param <- csaw::readParam(minq=30, forward = NULL, restrict = restrictChr)
	regionparam <- csaw::readParam(minq=30, restrict = restrictChr)

	# Count reads into sliding windows
	data <- csaw::strandedCounts(bam.files, param=param, ext=frag.len, width=win.width, bin = TRUE)
	colnames(data) <- rownames(design)
	SummarizedExperiment::colData(data) <- c(SummarizedExperiment::colData(data), design)

	# Get counts for 2kb local region surrounding each bin
	surrounds <- 2000
	neighbor <- suppressWarnings(
		GenomicRanges::trim(
			GenomicRanges::resize(SummarizedExperiment::rowRanges(data),
							   surrounds, fix = "center")
			))

	wider <- suppressWarnings(
		csaw::regionCounts(bam.files, param = regionparam, regions = neighbor, ext = frag.len)
	)

	colnames(wider) <- rownames(design)
	SummarizedExperiment::colData(wider) <- c(SummarizedExperiment::colData(wider), design)

	## take out groups --> Generate filter statistics for each group (based on local enrichment)
	filterstat <- lapply(unique(design$group), function(x){
		stat <- csaw::filterWindows(data[, data$group == x],
					    wider[, wider$group == x],
					    type = "local")
		return(stat)
	})

	# Require X-fold enrichment over local background to keep the window (similar to MACS)
	keep <- lapply(filterstat, function(x) {
		kp <- x$filter > log2(foldChange)
		return(kp)
	})

	filtered.data <- lapply(keep, function(keep) {
		return(data[keep,])
	})

	## merge nearby windows (within 50bp) to get broader TSS
	merged <- lapply(filtered.data, function(d) {
		return(csaw::mergeWindows(d, tol = 10L, ignore.strand = FALSE))
	})

	## write merged output for each group
	message("Writing output .bed files per group")
	mapply(function(bedfile, group) {
		rtracklayer::export.bed(object = bedfile$region, con = group)
	}, bedfile = merged, group = paste0(outfile_prefix, "_" , design$group, ".bed") )

	## write out the union of merges
	message("Writing merged .bed files")
	merged <- lapply(merged, function(x) return(x$region))
	mergedall <- base::Reduce(S4vectors::union, merged)
	rtracklayer::export.bed(mergedall,  con = paste(outfile_prefix, "merged.bed", sep = "_"))

	# return data
	output <- list(counts.windows = data, counts.background = wider, filter.stats = filterstat)
	return(output)
}
