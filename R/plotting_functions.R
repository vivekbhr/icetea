
#' Plot read statistics from the CapSet object
#'
#' @param CapSet The \code{\link{CapSet}} object
#' @param plotType The type of plot to make. Choose from "numbers" or "proportions". If "proportions"
#'                 is selected, the proportion of reads w.r.t total demultiplexed reads per sample
#'                 would be plotted
#' @param outfile Output file name. (filename extention would be used to determine type).
#'                If outfile not specified, the plot would be retured on the screen
#'
#' @return A ggplot object, or a file. Plot showing the number/proportion of reads in each category, per sample
#' @export
#'
#' @examples
#' \dontrun{
#' plot_readStats(CapSet = cs, plotType = "numbers")
#' }

plot_readStats <- function(CapSet, plotType = c("numbers", "proportions"), outfile = NULL) {
	si <- sampleInfo(CapSet)
	## get info on how many columns present
	sicols <- colnames(si)
	fields_toplot <- c("demult_reads", "num_mapped", "num_filtered", "num_intss")
	msg <- "Plotting following information :"
	fields <- sapply(fields_toplot, function(x) if(x %in% sicols) {
		return(x)
	})
	message(cat(msg, fields))
	## prepare df
	si_stats <- data.frame(sample = si$samples,
			       demutiplexed_reads = si$demult_reads,
			       mapped_reads = si$num_mapped,
			       duplicate_free_reads = si$num_filtered,
			       reads_within_TSS = si$num_intss)
	if(plotType == "proportions") {
		si_stats[-1] <- si_stats[-1]/si_stats$demutiplexed_reads
		y_label <- "Proportion of demultiplexed reads"
	} else {
		y_label <- "Number of reads"
	}

	si_stats <- reshape2::melt(si_stats, id.vars = "sample")
	# plot stacked barchart
	p <- ggplot(si_stats, aes_string("sample", "value", fill = "variable")) +
		geom_bar(stat = "identity", position = "stack") +
		theme_light(base_size = 16)  +
		scale_fill_brewer(type = "seq", palette = "Blues") +
		coord_flip() +
		labs(x = "Sample", y = y_label, fill = "Category")
	# return
	if(!(is.null(outFile))) {
		ggsave(outFile, plot = p, dpi = 300)
	} else {
		return(p)
	}

}



#' Plot precision of TSS calling between multiple samples
#'
#' @param TSSbedFiles BED files with TSS called by the CAGE wrapper
#' @param sampleNames Labels for Input samples
#' @param reference Reference Transcrips/Genes as a \code{\link{GRanges}} object
#' @param distanceCutoff Maximum distance (in base pairs) from reference TSS to plot
#' @param outFile Output file name (filename extention would be used to determine type).
#'                If outfile not specified, the plot would be retured on the screen
#'
#' @return A ggplot object, or a file. Plot showing % TSS called per sample w.r.t their distance to TSS
#' @export
#' @importFrom ggplot2 aes_string stat_ecdf theme_light scale_x_continuous scale_color_brewer ggsave
#'
#' @examples
#' \dontrun{
#' library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#' transcripts <- transcripts(dm6GTF)
#' files <- system.file("extdata", "testTSS.bed", package = "mapcapR")
#' plot_TSSprecision(TSSbedFiles = files, sampleNames = "testTSS",
#' 		reference = transcripts , distanceCutoff = 500,
#' 		outFile = "TSS_detection_precision.png")
#' }

plot_TSSprecision <- function(TSSbedFiles, sampleNames, reference, distanceCutoff = 500 , outFile = NULL) {
	# resize gene/transcript file to start
	refRanges <- GenomicRanges::resize(reference, width = 1, fix = "start")
	refRanges <- unique(refRanges)
	# read bed files
	tssData <- lapply(TSSbedFiles, rtracklayer::import.bed)
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
	p <- ggplot(tssdistances, aes_string("distances", col = "sample")) +
		stat_ecdf(geom = "step", size = 1) +
		theme_light(base_size = 14)  +
		scale_x_continuous(limits = c(0,distanceCutoff)) +
		scale_color_brewer(palette = "Set2") +
		labs(x = "Distances from nearby TSS (in bp)", y = "Cumulative Fraction",
		     title = "TSS precisions", col = "Category")
	if(!(is.null(outFile))) {
		ggsave(outFile, plot = p, dpi = 300)
	} else {
		return(p)
	}


}
