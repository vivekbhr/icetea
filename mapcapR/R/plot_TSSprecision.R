
#' Plot precision of TSS calling between multiple samples
#'
#' @param TSSbedFiles BED files with TSS called by the CAGE wrapper
#' @param sampleNames Labels for Input samples
#' @param reference Reference Transcrips/Genes GRanges
#' @param distanceCutoff Maximum distance upto which to plot precision
#' @param outFile Output png file name
#'
#' @return A plot showing % TSS called per sample w.r.t their distance to TSS
#' @export
#'
#' @examples
#' \dontrun{
#' library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#' transcripts <- transcripts(dm6GTF)
#' files <- system.file("extdata", c("tssclusters1.bed","tssclusters2.bed"), package = "FATSCapR")
#' plot_TSSprecision(TSSbedFiles = files, sampleNames = c("test1", "test2"),
#' 		reference = transcripts , distanceCutoff = 500,
#' 		outFile = "TSS_detection_precision.png")
#' }

plot_TSSprecision <- function(TSSbedFiles, sampleNames, reference, distanceCutoff = 500 , outFile = NULL) {
	# resize gene/transcript file to start
	refRanges <- GenomicRanges::resize(reference, width = 1, fix = "start")
	refRanges <- unique(refRanges)
	# read bed files
	tssData <- lapply(TSSbedFiles, bedToGRanges)
	names(tssData) <- sampleNames

	# get distances of bed entries from nearest TSS
	tssdistances <- sapply(tssData, function(x){
		y <- GenomicRanges::distanceToNearest(x,refRanges)
		return(as.data.frame(y)$distance)
	})

	# melt to df
	tssdistances <- plyr::ldply(tssdistances, data.frame)
	colnames(tssdistances) <- c("sample","distances")

	# plot ECDF with distance cutoff
	png(outFile, res = 100, width = 800, height = 800)
	ggplot(tssdistances, aes(distances, col = sample)) +
		stat_ecdf(geom = "step", size = 1) +
		theme_light(base_size = 14)  +
		scale_x_continuous(limits = c(0,distanceCutoff)) +
		scale_color_brewer(palette = "Set2") +
		labs(x = "Distances from nearby TSS (in bp)", y = "Cumulative Fraction",
		     title = "TSS precisions", col = "Category")
	dev.off()
}
