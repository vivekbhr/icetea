
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
#' ?plot_TSSprecision
#'

plot_TSSprecision <- function(TSSbedFiles, sampleNames, reference, distanceCutoff = 500 , outFile = NULL) {
	# resize gene/transcript file to start
	refRanges <- GenomicRanges::resize(reference, width = 1, fix = "start")
	# read bed files
	tssData <- lapply(TSSbedFiles, bedToGRanges)
	names(tssData) <- sampleNames

	# get distances of bed entries from nearest TSS
	tssdistances <- lapply(tssData, function(x){
		y <- GenomicRanges::distanceToNearest(x,refRanges)
		return(as.data.frame(y))
	})
	lapply(tssdistances, function(x) print(head(x)))
	# take cumulative sum of distance from nearest TSS
	tssdistances <- lapply(tssdistances, function(x){
		dist <- dplyr::arrange(x,distance)
		dist$numbers <- 1:nrow(dist)
		return(dist)
	})

	# plot line chart of num of detected TSS vs distance from TSS
	png(outFile, res = 100, width = 800, height = 800)
	lapply(1:length(tssdistances), function(n){
		# get name and color
		name <- names(tssdistances)[n]
		colpal <- colorRampPalette(RColorBrewer::brewer.pal(7,"Set1"))(length(tssdistances))
		test <- tssdistances[[n]]
		test[test$distance >  distanceCutoff, "distance"] <-  distanceCutoff
		test$numbers <- test$numbers*100/max(test$numbers) # convert to %age
		if(n == 1){
			plot(numbers ~ distance, test, type = "l", lwd = 2, col = colpal[n], xaxt = "n",
			     xlab = "Distance from annotated TSS", ylab= "% of known TSS detected",
			     main = "TSS detection precision")
			axis(side = 1, at = seq(0,distanceCutoff,50))
			legend(10,100, legend = name, fill = colpal[n], cex = 0.7)
		} else {
			lines(numbers ~ distance, test, type = "l",lwd = 2, col = colpal[n])
			legend(n*50,100, legend = name, fill = colpal[n], cex = 0.7)
		}
	})
	dev.off()
}
