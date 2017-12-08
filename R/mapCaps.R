
# check for trimmed barcodes present in the fastq header (RAMPAGE/MAPCap)
.checkTrimBarcodes <- function(fastq) {
	message("Checking for fastq file")
	stopifnot(file.exists(fastq))
	# test for trimmed R2 index
	message("Checking for trimmed barcodes")
	read <- gzfile(fastq)
	data <- readLines(read,100)
	close(read)

	# Check if the read header has "#" (which is introduced during trimming.)
	header <- data[seq(1,100,4)]
	if (unique(grepl("#",header)) != TRUE) {
		stop("Stop! fastq reads seem untrimmed. Run trimFastqIndex first.")
	}

}

#' Map the CAGE data
#'
#' @param CapSet An object of class \code{\link{CapSet}}
#' @param genomeIndex path to the Subread index file. Should end with the basename of the index.
#' @param outdir output directory path
#' @param nthreads number of threads to use for mapping.
#' @param logfile a log file to write the processing message.
#' @param ... additional arguments passed to the \code{\link{subjunc}} function.
#'
#' @return modified CapSet object with mapping information. Mapped and sorted BAM files are saved in `outdir`.
#' @export
#'
#' @examples
#'
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

	# test for trimmed R2 index
	if (expMethod %in% c("RAMPAGE", "MAPCap")) {
		.checkTrimBarcodes(CapSet@trimmed_R2)
	}

	if (demult) {
		samplelist <- sampleInfo$samples
		R1_list <- as.character(sampleInfo$demult_R1)
		R2_list <- as.character(sampleInfo$demult_R2)
	} else {
		samplelist <- list("trimmed")
		R1_list <- CapSet@trimmed_R1
		R2_list <- CapSet@trimmed_R2
	}

	# open a logfile if given
#	if(!is.null(logfile)){
#		log_con <- file(logfile)
#	} else log_con <- NULL

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

	# Close logfile
	#if (!is.null(logfile)) close(log_con)

	# edit sampleinfo of CapSet
	si <- sampleInfo(CapSet)
	maptable <- as.data.frame(t(mapstat[c(1,3),]))
	si$mapped_file <- as.character(maptable$Samples)
	si$num_mapped <- as.numeric(maptable$NumMapped)
	sampleInfo(CapSet) <- si

	return(CapSet)
}



