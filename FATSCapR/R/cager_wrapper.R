
#' A wrapper for CAGE Normalization and TSS calling
#'
#' @param input path to input bam file(s)
#' @param labels sample labels
#' @param ncores number of threads (leave NULL for one)
#' @param tpmCutoff TPM cutoff for TSS calling
#' @param genome The name of corresponding BSGenome object (must be installed)
#'
#' @return QC plots, Normalized bedgraph files and TSS bed files
#' @export
#'
#' @examples
#' 
#' cager_wrapper(input, labels, ncores = NULL, tpmCutoff = 10, genome = "BSgenome.Dmelanogaster.UCSC.dm6")
#' 

cage_wrapper <- function(input, labels, ncores = NULL, tpmCutoff = 10, genome = "BSgenome.Dmelanogaster.UCSC.dm6"){
	# build mycage
	mycage <- new("CAGEset", genomeName = genome,
			  inputFiles = input, inputFilesType = "bam", sampleLabels = labels)
	ctss <- getCTSS(mycage, removeFirstG = TRUE, correctSystematicG = TRUE)
	
	# plot correlation
	corr.m <- plotCorrelation(mycage, samples = "all", method = "pearson")
	# get lib sizes
	print("Library Sizes : ")
	print(librarySizes(mycage))
	# plot rev cums
	plotReverseCumulatives(mycage, fitInRange = c(5, 1000), onePlot = TRUE)
	
	# After loking into the plot I can determine suitable alpha and T for normalization
	# I can use it if i want to normalize by power law, but disabled by default
	#normalizeTagCount(mycage, method = "powerLaw", fitInRange = c(5, 50000), alpha = 1.05, T = 1*10^4)
	
	normalizeTagCount(mycage, method = "simpleTpm")
	# export normalized data to bedgraph
	exportCTSStoBedGraph(mycage, values = "normalized", oneFile = FALSE)
	
	# Check if multicore needed
	multi <- ifelse(is.null(ncores), FALSE, TRUE)
	
	if (multi == TRUE) {
		PARAM = BiocParallel::MulticoreParam(ncores)
	} else {
		PARAM = BiocParallel::MulticoreParam(1)
	}
	
	# TSS clustering
	clusterCTSS(object = mycage, threshold = tpmCutoff, thresholdIsTpm = TRUE, nrPassThreshold = 1, 
			method = "distclu", maxDist = 20, removeSingletons = TRUE, keepSingletonsAbove = 10,
			useMulticore = multi, nrCores = ncores)
	tc <- BiocParallel::bplapply(sampleLabels(mycage), function(x) tagClusters(mycage, sample = x), BPPARAM = PARAM)
	
	# get TSS width (combining all samples)
	message("Annotating TSS")
	cumulativeCTSSdistribution(mycage, clusters = "tagClusters")
	quantilePositions(mycage, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
	tc <- BiocParallel::bplapply(sampleLabels(mycage), function(x) {
		tagClusters(mycage, sample = x, returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
	}, BPPARAM = PARAM)
	
	# export bed files
	exportToBed(object = mycage, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneFile = FALSE, colorByExpressionProfile = TRUE)
	
}
