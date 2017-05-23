
#' Annotate detected Transcription Start Sites
#'
#' @param tssFile BED file with detected TSS/differential TSS results.
#' @param txdb A txdb object.
#' @param plot Type of plot to make (choose from "numbers", "proportions" or NA for no plot)
#'
#' @return Annotation of detected TSS
#' @export
#'
#' @examples
#'
#'

annotate_TSS <- function(tssFile, txdb, plot = NA) {

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

	## plot if asked
	if(!is.na(plot)) {
		if(plot == "numbers") {
			t <- as.data.frame(apply(x[-1], 1, function(x) (x/sum(x))*100 ) )
			colnames(t) <- x$.id
			t$Feature <- rownames(t)
			t <- melt(t)
			n <- "% "
		} else if(plot == "proportions") {
			t <- melt(x)
			n <- "Number "
		} else {
			warning("Plot type neither 'number' nor 'proportion'.")
		}

		print(ggplot(t, aes(Feature, value, fill = variable)) +
			geom_bar(stat = "identity", position = "dodge") +
			scale_fill_brewer(palette = "Set1") +
			labs(x = "Feature", y = paste0(n, "of TSS"), fill = "Sample")
		)
	}

	return(final)

}


## resolve 1:many mapping issue by prioritising some features over others
rankdf <- data.frame(feature = c("fiveUTR","promoter", "intron","coding","spliceSite","threeUTR","intergenic"),
		     rank = c(1,2,3,4,5,6,7))

#' Assign feature ranks on a VariantAnnotation output
#'
#' @param x output from VariantAnnotation
#'
#' @return A list of ranks
#'
#'
getranks <- function(x) {
	x$rank <- sapply(x$LOCATION, function(y) {
		return(rankdf[rankdf$feature == y, "rank"])
	})
	return(x)
}

#' Get features with the best rank for each TSS
#'
#' @param x output of getranks
#'
#' @return A data frame with counts
#'
#'
splitranks <- function(x) {
	l <- lapply(split(x, x$QUERYID), unique)
	l2 <- lapply(l, function(y) {
		return(y[which(y$rank == min(y$rank)),])
	})
	l3 <- plyr::ldply(l2, data.frame)
	return(l3)
}


#' Melt the output df from splitranks
#'
#' @param x
#'
#' @return
#'
#' @examples
#'
#'
melt <- function(x) {
	vars <- colnames(x[2:ncol(x)])
	d <- plyr::unrowname(reshape(x, direction = "long", idvar = ".id", varying = vars) )
	d$time <- vars[d$time]
	colnames(d) <- c("variable", "Feature", "value")
}

