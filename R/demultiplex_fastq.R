
#' Filter a read using the borcode in the header
#'
#' @param idx_name barcode ID (character)
#' @param fq_id fastq read header
#' @param maxM max allowd mismatches
#'
#' @return logical vector (which IDs to keep)
#'
filter_byIDx <- function(idx_name, fq_id, maxM){
	idx_name <- Biostrings::DNAString(idx_name)
	sep1 <- vapply(strsplit(as.vector(fq_id), "#"), "[[", character(1), 2)
	sep2 <- vapply(strsplit(sep1, ":"), "[[", character(1), 1)
	sep2 <- Biostrings::DNAStringSet(sep2)
	keep <- as.logical(Biostrings::vcountPattern(idx_name, sep2, max.mismatch = maxM))
	return(keep)
}


#' Split paired-end fastq by barcodes
#'
#' @param idx_name barcode ID
#' @param outfile_R1 output fastq file : Read 1
#' @param outfile_R2 output fastq file : Read 2
#' @param fastq_R1 input fastq file : Read 1
#' @param fastq_R2 input fastq file : Read 2
#' @param max_mismatch max allowd mismatches
#'
#' @return split fastq files
#'
split_fastq <- function(idx_name, outfile_R1, outfile_R2, fastq_R1, fastq_R2, max_mismatch) {

	## open input stream
	stream_R1 <- ShortRead::FastqStreamer(fastq_R1)
	stream_R2 <- ShortRead::FastqStreamer(fastq_R2)
	on.exit(close(stream_R1))
	on.exit(close(stream_R2), add = TRUE)

	kept_reads <- 0
	## filter fastq
	repeat {
		# input chunk
		fq_R1 <- ShortRead::yield(stream_R1)
		fq_R2 <- ShortRead::yield(stream_R2)
		if (length(fq_R1) == 0) {
			break
		}

		# get R1 and R2 reads and quality
		fq_R1read <- ShortRead::sread(fq_R1)
		fq_R1qual <- Biostrings::quality(fq_R1)
		fq_R2read <- ShortRead::sread(fq_R2)
		fq_R2qual <- Biostrings::quality(fq_R2)

		# trim barcodes from R2 and prepare header
		# copy the index sequence (position 6 to 12)
		sample_idx <- IRanges::narrow(fq_R2read, 6, 12)
		# copy pcr barcode (pos 1 to 5 + pos 12 to 13)
		umi_barcodes <- Biostrings::DNAStringSet(paste0(IRanges::narrow(fq_R2read, 1, 5),
								IRanges::narrow(fq_R2read, 12, 13) ) )
		#copy replicate demultiplexing barcode  (pos 3 to 4)
		rep_idx <- IRanges::narrow(fq_R2read, 3, 4)

		barcode_string <- paste(sample_idx, umi_barcodes, rep_idx, sep = ":")
		# now trim the first 13 bp off the read and quality
		fq_R2read <- IRanges::narrow(fq_R2read, 14, width(fq_R2) )
		fq_R2qual <- IRanges::narrow(fq_R2qual, 14, width(fq_R2) )

		## make a new ShortReadQ object with new header and trimmed reads
		# new fastq R2
		fq_R2new <- ShortRead::ShortReadQ(fq_R2read,
						  fq_R2qual,
						  Biostrings::BStringSet(paste(id(fq_R2), barcode_string, sep = "#")) )

		# new fastq R1 (seq/qual not modified, just barcodes copied from R2)
		fq_R1new <- ShortRead::ShortReadQ(fq_R1read,
						  fq_R1qual,
						  Biostrings::BStringSet(paste(id(fq_R1), barcode_string, sep = "#")) )

		# demultiplex using given sample barcodes
		idx_name <- Biostrings::DNAString(idx_name)
		sample_idx <- Biostrings::DNAStringSet(sample_idx)
		id2keep <- as.logical(Biostrings::vcountPattern(idx_name, sample_idx, max.mismatch = max_mismatch))

		#id2keep <- filter_byIDx(idx_name,
		#			fq_id = ShortRead::id(fq_R2),
		#			maxM = max_mismatch)
		# append to destination
		ShortRead::writeFastq(fq_R1new[id2keep], outfile_R1, "a")
		ShortRead::writeFastq(fq_R2new[id2keep], outfile_R2, "a")
		# add to count
		kept_reads <- kept_reads + sum(id2keep)
	}

	return(kept_reads)
}


#' Demultiplex (tagged) fastq files using sample barcodes
#'
#' @param CapSet CapSet object created using \code{\link{newCapSet}} function
#' @param max_mismatch maximum allowd mismatches
#' @param outdir path to output directory
#' @param nthreads No. of threads to use
#'
#' @return de-multiplxed fastq files corresponding to each barcode. The files are written
#'         on disk with the corresponding sample names as specified in sampleBarcodes(CapSet)
#' @export
#'
#' @examples
#' \dontrun{
#' idxlist <- c("CAAGTG", "TTAGCC", "GTGGAA", "GGTTAC", "TGTGAG", "CATCAC")
#' fnames <- c("GFPa","GFPb","BEAF32a","BEAF32b","BEAF32Rrp6a","BEAF32Rrp6b")
#'
#' cs <- newCapSet(expMethod = 'MAPCap', fastqType = 'paired',
#'		fastq_R1 = '../../2017_MAPCap/MAPCap_F/01_fastq/trimmed/MAPCap_F_R1.fastq.gz',
#'		fastq_R2 = '../../2017_MAPCap/MAPCap_F/01_fastq/trimmed/MAPCap_F_R2.fastq.gz',
#'		sampleInfo = data.frame(row.names = idxlist, files = fnames))
#'
#' system.time(demultiplex_fastq(cs, max_mismatch = 2))
#' }

demultiplex_fastq <- function(CapSet, max_mismatch, outdir, nthreads = 1) {

	sampleinfo <- sampleInfo(CapSet)
	destinations <- as.character(sampleinfo[,1])
	idx_list <- as.character(rownames(sampleinfo))

	## get the fastq to split (raise error if fastq untrimmed/not existing)
	if (is.null(CapSet@trimmed_R1) | is.null(CapSet@trimmed_R2)) {
		stop("FASTQ files with trimmed headers absent, can't demultiplex")
	} else {
		.checkTrimBarcodes(CapSet@trimmed_R1)
		fastq_R1 <- CapSet@trimmed_R1
		fastq_R2 <- CapSet@trimmed_R2
	}

	param = BiocParallel::MulticoreParam(workers = nthreads)
	message("de-multiplexing the FASTQ file")
	## filter and write
	info <- BiocParallel::bplapply(seq_along(destinations), function(i){
		split1 <- file.path(outdir, paste0(destinations[i],"_R1.fastq.gz"))
		split2 <- file.path(outdir, paste0(destinations[i],"_R2.fastq.gz"))

		kept <- split_fastq(idx_list[i], split1, split2,
				    fastq_R1, fastq_R2,
				    max_mismatch)
		dfout <- data.frame(R1 = split1, R2 = split2, kept_reads = kept)
		return(dfout)
		}, BPPARAM = param)

	## add post-demult info to sampleInfo
	keptinfo <- do.call(rbind, info)
	sampleinfo$demult_reads <- keptinfo$kept_reads
	sampleinfo$demult_R1 <- as.character(keptinfo$R1)
	sampleinfo$demult_R2 <- as.character(keptinfo$R2)

	## return object
	sampleInfo(CapSet) <- sampleinfo
	return(CapSet)

}
