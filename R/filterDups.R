
#' Filter PCR-duplicates from BAM file using internal UMIs
#'
#'
#' @param bamFile Input BAM file
#' @param outFile Output (filtered) BAM file
#'
#' @return Filtered BAM file (with only R1), after PCR duplicate removal
#' @export
#'
#' @examples
#'
#' bam <- system.file("extdata", "test_mapped.bam", package = "mapcapR")
#' filterDuplicates(bamFile = bam, outFile = "test_rmDup.bam")
#'


filterDups <- function(bamFile, outFile) {
	message(paste0("Removing PCR duplicates : ", bamFile))

	sparam <- Rsamtools::ScanBamParam(what = c("qname", "rname", "pos", "isize", "qwidth", "mapq"),
					  flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
					  			      isFirstMateRead = TRUE) )
	bamdf <- Rsamtools::scanBam(bamFile, param = sparam)
	bamdf <- do.call(as.data.frame, bamdf)
	bamdf$qname <- as.character(bamdf$qname)

	## assuming a bamdf available, this is what a filterfunc looks like
	filterDups <- function(bamdf) {
		bamdf$qname <- as.character(bamdf$qname)
		## split the df by chromosome into multiple df
		chroms <- factor(bamdf$rname, levels = unique(as.character(bamdf$rname)))
		bam2 <- S4Vectors::split(bamdf, chroms)

		## get duplicate stats for a given df
		getdupstats <- function(bamdf) {
			## split the df by pos (get one list per pos)
			bamdf.bypos <- split(bamdf, list(bamdf$pos))
			## extract umis from each df
			getumi <- function(x) {
				hdr <- vapply(strsplit(x$qname, "#"), "[[", character(1), 2)
				umi <- vapply(strsplit(hdr, ":"), "[[", character(1), 2)
				return(umi)
			}
			#getfraglength <- function(x) {
			#	fraglen <- (2*x$qwidth)  x$isize
			#	return(fraglen)
			#}
			umis <- lapply(bamdf.bypos, getumi)
			#fraglengths <- lapply(bamdf.bypos, getfraglength)

			dupStats_umi <- unlist(lapply(umis, function(x) !(duplicated(x)) ))
			#dupStats_fraglen <- unlist(lapply(fraglengths, function(x) !(duplicated(x)) ))

			# final dupstats (both UMI and fragment length are same --> remove reads, else keep)
			#dupStats <- dupStats_umi | dupStats_fraglen
			return(dupStats_umi)
		}

		## run for all
		dupstats_perchr <- lapply(bam2, getdupstats)
		return(unlist(dupstats_perchr))

	}

	## create rule
	rule <- S4Vectors::FilterRules(list(filterDups))

	## filter command
	Rsamtools::filterBam(file = bamFile,
			     destination = outFile,
			     filter = rule, param = sparam )

}


#' Filter PCR-duplicates from mapped files using internal UMIs
#'
#' @description This script considers the read mapping start position and the UMI to determine whether a
#'              read is a PCR duplicate. All PCR duplicates are then removed and one entry per read is kept.
#'              In case of paired-end reads (MAPCap/RAMPAGE), only one end (R1) is kept after filtering.
#' @param CapSet an object of class \code{\link{CapSet}}
#' @param outdir output directory for filtered BAM files
#'
#' @return modified CapSet object with filtering information. Filtered BAM files are saved in `outdir`.
#' @export
#'
#' @examples
filterDuplicates <- function(CapSet, outdir) {
	si <- sampleInfo(CapSet)
	bamfiles <- si$mapped_file

	# first check if the bam files exist
	lapply(bamfiles, function(f) {
		if (!(file.exists(f) )) stop(paste0("mapped file ", f, " doesn't exist!",
						  "Please update your Capset object with valid file paths ",
						  "using sampleInfo(CapSet). "))
	})
	# then prepare outfile list
	outfiles <- file.path(outdir, paste0(si$samples, ".filtered.bam") )
	# run the filter duplicates function on all files
	mapply(filterDups, bamfiles, outfiles)

	# collect post-filtering stats
	mapstat <- sapply(outfiles, Rsubread::propmapped, simplify = TRUE)
	# update CapSet
	maptable <- as.data.frame(t(mapstat[c(1,3),]))
	si$filtered_file <- as.character(maptable$Samples)
	si$num_filtered <- as.numeric(maptable$NumMapped)
	sampleInfo(CapSet) <- si

	validObject(CapSet)
	return(CapSet)

}
