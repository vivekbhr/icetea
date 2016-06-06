##
##	GENERIC FUNCTIONS FOR INTERNAL USE IN FATSCapR package
##


#' Import TSS bed file as GRanges
#'
#' @param bed BedFile to import
#'
#' @return GRanges object
#' @export
#'
#' @examples
#' bedToGRanges(bed)
#'

bedToGRanges <- function(bed){
	sbed <- read.delim(bed,skip=1, header = FALSE)
	sbed <- GRanges(seqnames = sbed$V1, ranges = IRanges(start = sbed$V2, end = sbed$V3), strand = sbed$V6)
	return(sbed)
}
