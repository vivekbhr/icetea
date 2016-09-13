
#' A wrapper for CAGE Normalization and TSS calling
#'
#' @param input path to input bam file(s)
#' @param labels sample labels
#' @param ncores number of threads (leave NULL for one)
#' @param tpmCutoff TPM cutoff for TSS calling
#' @param genome The name of corresponding BSGenome object (must be installed)
#' @param find_TSSshift logical, indicating if TSS shift between two samples are to be calculated
#' @param scoreshift_groupX Test sample for TSS shift
#' @param scoreshift_groupY Control sample for TSS shift
#' @param promoterShift_outFile Output file name for calculated TSS shift positions
#'
#' @return QC plots, Normalized bedgraph files and TSS bed files
#' @export
#'
#' @examples
#'
#' bam <- system.file("extdata", "test_mapped.bam", package = "FATSCapR")
#' cager_wrapper(input = bam, labels = "test", tpmCutoff = 10, genome = "BSgenome.Dmelanogaster.UCSC.dm6")
#'

cage_wrapper <- function(input, labels, ncores = NULL, tpmCutoff = 10, genome = "BSgenome.Dmelanogaster.UCSC.dm6",
			 find_TSSshift = FALSE, scoreshift_groupX, scoreshift_groupY, promoterShift_outFile){
	## 1... build mycage
	mycage <- new("CAGEset", genomeName = genome,
			  inputFiles = input, inputFilesType = "bam", sampleLabels = labels)
	tryCatch({
		ctss <- CAGEr::getCTSS(mycage, removeFirstG = TRUE, correctSystematicG = TRUE)
	}, error = function(e){
		if(grepl("negative",e) == TRUE){
			print("It seems that the file is too big. Try splitting the file by chromosome and
				process one by one.")
		} else {
			print("Error: ", e)
		}
	})

	# plot correlation
	corr.m <- CAGEr::plotCorrelation(mycage, samples = "all", method = "pearson")
	# get lib sizes
	message("Library Sizes : ")
	libsize <- CAGEr::librarySizes(mycage)
	print(as.data.frame(libsize))

	# plot rev cums(only to check if sample follows power law)
	CAGEr::plotReverseCumulatives(mycage, fitInRange = c(5, 1000), onePlot = TRUE)
	CAGEr::normalizeTagCount(mycage, method = "simpleTpm")
	# export normalized data to bedgraph
	CAGEr::exportCTSStoBedGraph(mycage, values = "normalized", oneFile = FALSE)

	## 2 a. Prepare for step 2 and 3

	# Check if multicore needed
	multi <- ifelse(is.null(ncores), FALSE, TRUE)

	if (multi == TRUE) {
		PARAM = BiocParallel::MulticoreParam(ncores)
	} else {
		PARAM = BiocParallel::MulticoreParam(1)
	}
	## get some default thresholds
	message("Min quantile, Max Quantile, Max distance for clustering, sigletons to exclude, FDR cutoff")
	(qmin <- 0.1)
	(qmax <- 0.9)
	(clusterMaxdist <- 20)
	(sigletonCutoff <- 10)
	(fdr <- 0.01)
	shift <- all(!(is.null(scoreshift_groupX)),
		     !(is.null(scoreshift_groupY)),
		     !(is.null(promoterShift_outFile)) )

	## 2b. .. TSS clustering
	message("Clustering tags and annotating TSS")
	tryCatch({
		CAGEr::clusterCTSS(object = mycage, threshold = tpmCutoff, thresholdIsTpm = TRUE, nrPassThreshold = 1,
					 method = "distclu", maxDist = clusterMaxdist,
					 removeSingletons = TRUE, keepSingletonsAbove = sigletonCutoff,
					 useMulticore = multi, nrCores = ncores)
		tc <- BiocParallel::bplapply(sampleLabels(mycage), function(x) tagClusters(mycage, sample = x),
						     BPPARAM = PARAM)

		# get TSS width (combining all samples)
		message("Annotating TSS")
		CAGEr::cumulativeCTSSdistribution(mycage, clusters = "tagClusters")
		CAGEr::quantilePositions(mycage, clusters = "tagClusters", qLow = qmin, qUp = qmax)
		tc <- BiocParallel::bplapply(sampleLabels(mycage), function(x) {
			tagClusters(mycage, sample = x, returnInterquantileWidth = TRUE, qLow = qmin, qUp = 0.9)
		}, BPPARAM = PARAM)

		# export bed files
		CAGEr::exportToBed(object = mycage, what = "tagClusters", qLow = qmin, qUp = qmax,
					 oneFile = FALSE, colorByExpressionProfile = TRUE)

	},  error = function(e){
		print("Error: ", e)
		print("CAGEset object with normalized counts saved in current directory as : myCAGEset.Rdata")
		save(ctss, file= "myCAGEset.Rdata")
	})

	## 3..Aggregate TSS by positions -- > cluster --> promoter shifts

	tryCatch({
		# consensus clusters
		message(paste0("Aggregating TSS by positions. Max distance = ",clusterMaxdist*5) )
		CAGEr::aggregateTagClusters(mycage, tpmThreshold = tpmCutoff, qLow = qmin, qUp = qmax, maxDist = clusterMaxdist*5)
		CAGEr::quantilePositions(mycage, clusters = "consensusClusters", qLow = qmin, qUp = qmax)

		consensusCl <- lapply(labels, function(sample) {
			return(CAGEr::consensusClusters(mycage, sample = sample,
						 returnInterquantileWidth = TRUE, qLow = qmin, qUp = qmax) )
		})
		print("Writing out consensus clusters..")
		names(consensusCl) <- labels
		lapply(labels, function(x) {
			write.table(consensusCl[[x]], file = paste0(x,"_consensusClusters.tsv"), sep = "\t", row.names = FALSE)
		})

		## expression profiles
		message("Clustering TSS by expression.")
		CAGEr::getExpressionProfiles(mycage, what = "consensusClusters", tpmThreshold = tpmCutoff,
				      nrPassThreshold = 1, method = "som", xDim = 4, yDim = 2)
		CAGEr::plotExpressionProfiles(mycage, what = "consensusClusters")
		CAGEr::exportToBed(mycage, what = "consensusClusters", colorByExpressionProfile = TRUE)

		## shifting promoters
		if(isTRUE(shift)) {
			message(paste0("Calculating promoter shifts. Sample : ", scoreshift_groupX, " over : ", scoreshift_groupY) )

			CAGEr::cumulativeCTSSdistribution(mycage, clusters = "consensusClusters")
			CAGEr::scoreShift(mycage, groupX = scoreshift_groupX, groupY = scoreshift_groupY, testKS = TRUE, useTpmKS = FALSE)
			shifting.promoters <- CAGEr::getShiftingPromoters(mycage,tpmThreshold = tpmCutoff,
								   scoreThreshold = 0.6, fdrThreshold = fdr)
			write.table(shifting.promoters$total, promoterShift_outFile, sep = "\t", row.names = FALSE)
		} else {
			print("Skipping calculation of promoter shifts.")
		}
	})

}
