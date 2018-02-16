
## Filter duplicates from a data.frame, using position and UMI seq
filterdups_func <- function(bamdf) {
	bamdf$qname <- as.character(bamdf$qname)
	bamdf$pos <- as.numeric(bamdf$pos)
	## split the df by chromosome into multiple df
	chroms <- factor(bamdf$rname, levels = unique(as.character(bamdf$rname)))
	bam2 <- S4Vectors::split(bamdf, chroms)

	## get duplicate stats for a given df
	getdupstats <- function(bamdf) {
		## split the df by pos (get one list per pos)
		pos <- factor(bamdf$pos, levels = unique(as.numeric(bamdf$pos)))
		bamdf.bypos <- S4Vectors::split(bamdf, pos)
		## extract umis from each df
		getumi <- function(x) {
			hdr <- vapply(strsplit(x$qname, "#"), "[[", character(1), 2)
			umi <- vapply(strsplit(hdr, ":"), "[[", character(1), 2)
			return(!(duplicated(umi)))
		}
		dupStats_umi <- lapply(bamdf.bypos, getumi)
		return(dupStats_umi)
	}
	#getfraglength <- function(x) {
	#	fraglen <- (2*x$qwidth)  x$isize
	#	return(fraglen)
	#}
	#dupStats_fraglen <- unlist(lapply(fraglengths, function(x) !(duplicated(x)) ))
	# final dupstats (both UMI and fragment length are same --> remove reads, else keep)
	#dupStats <- dupStats_umi | dupStats_fraglen

	## run for all
	dupstats_perchr <- lapply(bam2, getdupstats)
	return(unlist(dupstats_perchr))

}


#' Filter PCR-duplicates from BAM file using internal UMIs
#'
#'
#' @param bamFile Input BAM file
#' @param outFile Output (filtered) BAM file
#'
#' @return Filtered BAM file (with only R1), after PCR duplicate removal
#'

filterDups <- function(bamFile, outFile) {
	message(paste0("Removing PCR duplicates : ", bamFile))
	# get baminfo
	sparam <- Rsamtools::ScanBamParam(what = c("qname", "rname", "pos"), #, "isize", "qwidth", "mapq"
					  flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
					  			      isFirstMateRead = TRUE) )

	## create rule
	rule <- S4Vectors::FilterRules(list(filterdups_func))

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
#' @importFrom methods validObject
#' @export
#'
#' @examples
#'
#' # before running this
#' 1. Create a CapSet object
#' 2. de-multiplex the fastqs
#' 3. map them
#'
#' # load a previously saved CapSet object
#' dir <- system.file("extdata", package="icetea")
#' cs <- load(file.path(dir, "CSobject.Rdata"))
#'
#' # filter duplicate reads from mapped BAM files
#' cs <- filterDuplicates(cs, outdir = dir)
#'

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
