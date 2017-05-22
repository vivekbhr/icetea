
#' Annotate detected Transcription Start Sites
#'
#' @param tssFile BED file with detected TSS/differential TSS results.
#' @param txdb A txdb object.
#'
#' @return Annotation of detected TSS
#' @export
#'
#' @examples
#'
#'

annotate_TSS <- function(tssFile, txdb){

	# import TSS file
	x <- rtracklayer::import.bed(tssFile)
	# Annotate
	db <- VariantAnnotation::locateVariants(query = x,
						subject = txdb,
						VariantAnnotation::AllVariants(
							promoter =VariantAnnotation::PromoterVariants(upstream = 500,
												      downstream = 0)))
	## resolve 1:many mapping isues using ranks from rankdf
	t <- data.frame(QUERYID = db$QUERYID, LOCATION = db$LOCATION)
	tt <- getranks(t)
	ttt <- splitranks(tt)

	## Return a table of tss counts per feature
	final <- as.data.frame(table(ttt$LOCATION))
	return(final)

}


## resolve 1:many mapping issue by prioritising UTR over others
rankdf <- data.frame(feature = c("fiveUTR","promoter", "intron","coding","spliceSite","threeUTR","intergenic"),
		     rank = c(1,2,3,4,5,6,7))

getranks <- function(x) {
	x$rank <- sapply(x$LOCATION, function(y) {
		return(rankdf[rankdf$feature == y, "rank"])
	})
	return(x)
}

splitranks <- function(x) {
	l <- lapply(split(x, x$QUERYID), unique)
	l2 <- lapply(l, function(y) {
		return(y[which(y$rank == min(y$rank)),])
	})
	l3 <- plyr::ldply(l2, data.frame)
	return(l3)
}
