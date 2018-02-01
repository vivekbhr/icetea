
#' Detection of Trancription start sites based on local enrichment
#'
#' @param CapSet CapSet object created using \code{\link{newCapSet}} function
#' @param groups a character vector that contains group name of the sample, for replicate-based TSS
#'               calling (see example)
#' @param outfile_prefix Output name prefix for the .bed files
#' @param foldChange A fold change cutoff of local enrichment to detect the TSS. For samples with
#' 		'usual' amount of starting material and squencing depth (>=5ug starting material,
#' 		>= 5 mil reads/sample), a cut-off of 6 fold can be used. For samples with low
#' 		amount of material or sequencing depth, use a lower cut-off (eg. use 2-fold for
#' 		samples with 500ng starting material).
#' @param restrictChr Chromosomes to restrict the analysis to.
#'
#'
#' @return .bed files containing TSS position for each group, along with a bed file for consensus
#' 	   (union) TSS sites of all samples.
#' @export
#' @importFrom utils write.table
#'
#' @examples
#' \dontrun{
#' cs <- detect_TSS(Capset = cs, group = c("wt","wt","mut","mut"), outfile_prefix = "testTSS",
#'            foldChange = 6, restrictChr = c("2L","2R","X"))
#'}
#'


detect_TSS <- function(CapSet, groups,  outfile_prefix,
		       foldChange = 2, restrictChr = NULL) {

	# check whether group and outfile_prefix is provided
	if (missing(outfile_prefix)) stop("Please provide outfile_prefix!")
	if (missing(groups)) stop("Please provide groups!")

	# convert group to char
	si <- sampleInfo(CapSet)
	design <- data.frame(row.names = si$samples, group = as.character(groups) )

	if (is.null(si$filtered_file) ) {
		message("Filtered files not found under sampleInfo(CapSet). Using mapped files")
		bam.files <- si$mapped_file
	} else {
		bam.files <- si$filtered_file
	}
	if (sum(file.exists(bam.files)) != length(bam.files)) {
		stop("One or more bam files don't exist! Check sampleInfo(CapSet) ")
	}

	# Define read params
	frag.len <- NA
	win.width <- 10
	param <- csaw::readParam(minq = 2, forward = NULL, restrict = restrictChr)
	regionparam <- csaw::readParam(minq = 2, restrict = restrictChr)

	# Count reads into sliding windows
	data <- csaw::strandedCounts(bam.files, param = param, ext = frag.len, width = win.width, bin = TRUE)
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

	## merge nearby windows (within 10bp) to get broader TSS
	merged <- lapply(filtered.data, function(d) {
		return(csaw::mergeWindows(d, tol = 10L, ignore.strand = FALSE))
	})

	## write merged output for each group
	message("Writing output .bed files per group")
	mapply(function(bedfile, group) {
		rtracklayer::export.bed(object = bedfile$region, con = group)
	}, bedfile = merged, group = paste0(outfile_prefix, "_" , unique(design$group), ".bed") )

	## write out the union of merges
	message("Writing merged .bed files")
	merged <- lapply(merged, function(x) return(x$region))
	mergedall <- base::Reduce(S4Vectors::union, merged)
	rtracklayer::export.bed(mergedall,  con = paste(outfile_prefix, "merged.bed", sep = "_"))

	## Calculate prop reads in TSS per group
	message("Counting reads within detected TSS")
	reads_intss <- propReadsInBed(mergedall, bam.files)

	# Add the results as slots in CapSet
	Capset@counts.windows <- data
	Capset@counts.background <- wider
	CapSet@filter.stats <- S4Vectors::DataFrame(filterstat[[1]])

	#output <- list(counts.windows = data, counts.background = wider, filter.stats = filterstat)
	return(CapSet)
}

propReadsInBed <- function(regions, bams = NA) {
	counts <- GenomicAlignments::summarizeOverlaps(GenomicRanges::GRangesList(regions),
					     reads = Rsamtools::BamFileList(as.character(bams)),
					     mode = "Union",
					     inter.feature = FALSE)
	numreads <- SummarizedExperiment::assay(counts)
	return(t(numreads))
}
