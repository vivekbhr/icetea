
#' Detect differentially expressed Transcription Start Sites between two conditions (fit model)
#'
#' @param bam.files List of bam files to use
#' @param TSSfile A .bed file with TS sites to test. Normally it would be the *merged.bed file
#' 		  (output of \code{\link{detect_TSS}} command)
#'
#' @param design A data frame with rownames = sample names and a column called 'group'
#' 		that contains information about the sample group (see example)
#' @param outplots Output pdf filename for plots
#'
#' @param plotref Name of reference sample to plot for detection of composition bias in the
#' 		  data. Data is normalized using the TMM method to avoid composition bias.
#' @param testGroup Name of Test group (should be present in the design data frame)
#' @param contGroup Name of Control group (should be present in the design data frame)
#'
#' @return Returns an object of class DGEGLM.
#'
#' @export
#'
#' @examples
#'
#'

fit_diffTSS <- function(bam.files, TSSfile, design, outplots, plotref) {

	# first check if design df and bam files are accurate
	if(length(bam.files) != nrow(design)) {
		stop("Number of rows in design data frame doesn't match the number of bam files")
	} else {
		c <- design$sample
	}

	## Normalize for composition bias : TMM
	# important to try different bin sizes and see if the values are close to unity (low composition effect)
	restrictChr <- NULL
	regionparam <- csaw::readParam(minq=30, restrict = restrictChr)
	binned <- csaw::windowCounts(bam.files, bin = TRUE, width = 10000, param = regionparam)
	normfacs <- csaw::normOffsets(binned) # close to unity
	names(normfacs) <- c

	## visualize Effect of TMM normalization on composition bias
	y.bin <- csaw::asDGEList(binned)
	bin.ab <- edgeR::aveLogCPM(y.bin)
	adjc <- edgeR::cpm(y.bin, log=TRUE)
	colnames(adjc) <- c

	# plot ref sample vs all other samples
	message("plotting the composition effect")
	sampnumber <- ncol(adjc) - 1
	cols_toplot <- c[!(grepl(plotref, c))]
	n <- ceiling(sampnumber/3) #roundup to make divisible by 3

	par(cex.lab=1.5, mfrow=c(n,3))
	lapply(cols_toplot, function(x) {
		smoothScatter(bin.ab, adjc[,plotref] - adjc[,x], ylim = c(-6, 6),
			      xlab = "Average abundance", ylab = paste0("Log-ratio (", plotref," vs ", x,")") )

		abline(h = log2(normfacs[plotref]/normfacs[x]), col = "red")

	})

	# Import tss locations to test
	mergedall <- rtracklayer::import.bed(TSSfile)

	## get 5' read counts on the locations from the bam.files
	# function to resize reads
	ResizeReads <- function(reads, width=1, fix="start", ...) {
		reads <- as(reads, "GRanges")
		stopifnot(all(strand(reads) != "*"))
		resize(reads, width=width, fix=fix, ...)
	}

	# now read the data
	tsscounts <- GenomicAlignments::summarizeOverlaps(features = mergedall,
							  reads = bam.files,
							  preprocess.reads = ResizeReads)
	#### ------ Now do EdgeR ------ ####
	y <- csaw::asDGEList(tsscounts, norm.factors = normfacs)
	designm <- model.matrix(~0+group, design)
	y <- edgeR::estimateDisp(y, designm)
	fit <- edgeR::glmQLFit(y, designm, robust = TRUE)

	# check prior degrees of freedom (avoid Infs)
	message("Prior degrees of freedom : ")
	print(summary(fit$df.prior))

	pdf(outplots)

	## check that the Fit is good
	par(mfrow=c(1,2))
	# BCV plot
	o <- order(y$AveLogCPM)
	plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type = "l", lwd = 2,
	     ylim = c(0, 1), xlab = expression("Ave."~ Log[2] ~ "CPM"),
	     ylab = ("Biological coefficient of variation"))
	# dispersion plot
	edgeR::plotQLDisp(fit)

	## check with MDSplot if there is batch effect
	par(mfrow = c(2,2), mar = c(5,4,2,2))
	for (top in c(100, 500, 1000, 5000)) {
		out <- limma::plotMDS(edgeR::cpm(y, log = TRUE), main = top,
				      labels = design$group, top = top)
	}

	dev.off()

	## return the fit
	return(fit)

}



#' Detect differentially expressed Transcription Start Sites between two conditions (test)
#'
#' @param fit DGEGLM object (output of \code{\link{fit_diffTSS}} command )
#' @param testGroup Test group name
#' @param contGroup Control group name
#' @param tssFile The TSS .bed file used for \code{\link{fit_diffTSS}} command
#' @param MAplot_fdr FDR threshold to mark differentially expressed TSS in MAplot (NA = Don't make an MAplot)
#'
#' @return A \code{\link{GRanges}} object containing p-values of differential expression for each TSS.
#' @export
#'
#' @examples
#'
#'

detect_diffTSS <- function(fit, testGroup, contGroup, TSSfile, MAplot_fdr = NA) {

	# Import tss locations to test
	mergedall <- rtracklayer::import.bed(TSSfile)

	## Testing the differential TSS
	# make contrast matrix
	contr <- paste0("group", testGroup , "-group", contGroup)
	contrast <- limma::makeContrasts(contrasts = contr, levels = fit$design)
	# test
	results <- edgeR::glmQLFTest(fit, contrast = contrast)
	top <- as.data.frame(edgeR::topTags(results, n = Inf))

	# sort output by pvalue and return
	difftss <- mergedall[rownames(top) %>% as.numeric()]
	difftss$score <- top$FDR
	difftss$logFC <- top$logFC
	difftss$logCPM <- top$logCPM

	# MA plot
	if(!(is.na(MAplot_fdr)) ) {
	p <- ggplot(top, aes(logCPM, logFC, col = factor(top$FDR < MAplot_fdr))) +
			geom_point(alpha = 0.5) +
			geom_abline(slope = 0, intercept = 0) +
			scale_color_manual(values = c("grey40", "darkred")) +
			labs(col = "Differentially Expressed") +
			theme_gray(base_size = 14) +
			theme(legend.position = "top")
	print(p)
	}

	return(difftss)
}

