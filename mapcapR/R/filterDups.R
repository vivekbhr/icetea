
#' Filter PCR-duplicates from BAM file using internal UMIs
#'
#' @param bamFile Input BAM file
#' @param outFile Output BAM file
#'
#' @return Filtered BAM file (with only R1), after PCR duplicate removal
#' @export
#'
#' @examples
#'
#' bam <- system.file("extdata", "test_mapped.bam", package = "mapcapR")
#' filterDuplicates(bamFile = bam, outFile = "test_rmDup.bam")
#'


filterDuplicates <- function(bamFile, outFile) {

	sparam <- Rsamtools::ScanBamParam(what = c("qname", "flag", "rname", "pos", "mapq"),
			       flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
			       			      isFirstMateRead = TRUE) )
	# I also want hasUnmappedMate = TRUE eventually

	filterDups <- function(bam){

		# function to get umis from qname and bins
		dupumi_perbin <- function(qname, bins) {
			# Get the UMI sequence out from header
			getumi <- function(x) {
				hdr <- vapply(strsplit(x, "#"), "[[", character(1), 2)
				umi <- vapply(strsplit(hdr, ":"), "[[", character(1), 2)
				return(umi)
			}
			# Do it per bin
			umis <- lapply(split(qname, bins), getumi)
			# Get the duplicated umis per bin -> merge
			dupStatus <- lapply(umis, function(x) !(duplicated(x)) )
			return(unlist(dupStatus))
		}

		# split list by chromosome
		chroms <- factor(bam$rname, levels = unique(as.character(bam$rname)))
		bam2 <- S4Vectors::split(bam, chroms)

		getdupstats <- function(x){
			# round to nearest 1000th
			getEnd <- function(b) {
				maxp <- max(b$pos)
				y <- round(maxp, -3)
				z <- ifelse(maxp < y, y, y + 1000)
				return(z)
			}
			getStart <- function(b) {
				minp <- min(b$pos)
				y <- round(minp, -3)
				z <- ifelse(minp > y, y, y - 1000)
				return(z)
			}
			chromStart <- getStart(x)
			chromEnd <- getEnd(x)

			# make chrom-wise bins
			bins <- .bincode(x$pos, seq(chromStart,chromEnd, 1000))
			dupstats <- dupumi_perbin(as.character(x$qname), bins)
			return(dupstats)
		}

		dupstats <- lapply(bam2, getdupstats)
		return(unlist(dupstats))
	}

	rule <- S4Vectors::FilterRules(list(filterDups))

	# filter command
	Rsamtools::filterBam(file = bamFile, destination = outFile, filter = rule, param = sparam )

	return("Done!")
}
