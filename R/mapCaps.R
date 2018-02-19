#' Map the data from 5' profiling techniques
#'
#' @param CapSet An object of class \code{\link{CapSet}}
#' @param genomeIndex path to the Subread index file. Should end with the basename of the index.
#' @param outdir output directory path
#' @param nthreads number of threads to use for mapping.
#' @param logfile a log file to write the processing message.
#' @param ... additional arguments passed to the \code{\link{subjunc}} function.
#'
#' @return modified CapSet object with mapping information. Mapped and sorted BAM files are saved in `outdir`.
#'
#' @importFrom utils capture.output
#' @export
#'
#' @examples
#' \dontrun{
#' # before mapping :
#' # 1. Create a CapSet object
#' # 2. de-multiplex the fastqs
#'
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#'
#' # map the data (not available on windows)
#' library(Rsubread)
#' buildindex(basename = "dm6", reference = "/path/to/dm6genome.fa")
#' cs <- mapCaps(cs, genomeIndex = "dm6", outdir = dir, nthreads = 10)
#'
#' }
#'

mapCaps <- function(CapSet, genomeIndex, outdir, nthreads, logfile = NULL, ...){

	## extract info
	sampleInfo <- sampleInfo(CapSet)
	expMethod <- CapSet@expMethod

	# test whether the CapSet object has demultiplexed fastqs
	demult <- !(is.null(sampleInfo$demult_R1)) # TRUE/FALSE

	if (demult) {
		if (sum(sapply(sampleInfo$demult_R1, file.exists, simplify = TRUE)) != nrow(sampleInfo) ) {
			stop("One or more demultiplxed fastq files don't exist.
			     Please check file paths under sampleInfo(CapSet)")
		}
	}

	if (demult) {
		samplelist <- sampleInfo$samples
		R1_list <- as.character(sampleInfo$demult_R1)
		R2_list <- as.character(sampleInfo$demult_R2)
	} else {
		samplelist <- list("raw")
		R1_list <- CapSet@fastq_R1
		R2_list <- CapSet@fastq_R2
	}

	mapstat <- mapply(function(sample, R1, R2) {

		message(paste0("Mapping sample : ", sample))
		# Align using RSubread
		tmpout <- file.path(outdir, paste0(sample, ".tmp.bam") )

		capture.output(Rsubread::subjunc(index = genomeIndex,
				  readfile1 = R1,
				  readfile2 = R2,
				  output_file = tmpout,
				  nthreads = nthreads,
				  minFragLength = 10,
				  reportAllJunctions = TRUE,
				  ...),
		    file = logfile, append = TRUE)

		# Sort and Index
		message("Sorting and Indexing")
		dest <- file.path(outdir, sample)
		Rsamtools::sortBam(file = tmpout, destination = dest) # adds .bam suffix
		Rsamtools::indexBam(paste0(dest, ".bam"))
		file.remove(tmpout)

		# Get mapping stats
		stat <- Rsubread::propmapped(paste0(dest, ".bam"))
		return(stat)

	}, samplelist, R1_list, R2_list)

	# edit sampleinfo of CapSet
	si <- sampleInfo(CapSet)
	maptable <- as.data.frame(t(mapstat[c(1,3),]))
	si$mapped_file <- as.character(maptable$Samples)
	si$num_mapped <- as.numeric(maptable$NumMapped)
	sampleInfo(CapSet) <- si

	return(CapSet)
}
