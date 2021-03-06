#' Get data to create new ShortReadQ object after barcode trimming
#'
#' @param type character. expType of the CapSet object
#' @param fq_R1 character. fastq Read 1
#' @param fq_R2 character. fastq Read 2
#'
#' @importFrom ShortRead sread
#' @importFrom Biostrings quality DNAStringSet
#' @importFrom BiocGenerics width
#' @importFrom IRanges narrow
#'
#' @return A list with new R1 and R2 sequence, quality, barcode string and sample id
#'
get_newfastq <- function(type, fq_R1, fq_R2) {
    # get R1 and R2 reads and quality
    fq_R1read <- sread(fq_R1)
    fq_R1qual <- quality(fq_R1)
    fq_R2read <- sread(fq_R2)
    fq_R2qual <- quality(fq_R2)
    if (type == "MAPCap") {
        # trim barcodes (pos 1 to 13) from R2 and prepare header
        # copy the index sequence (position 6 to 11)
        sample_idx <- narrow(fq_R2read, 6, 11)
        # copy pcr barcode (pos 1 to 5 + pos 12 to 13)
        umi_barcodes <-
            DNAStringSet(paste0(
                narrow(fq_R2read, 1, 5),
                narrow(fq_R2read, 12, 13)
            ))
        #copy replicate demultiplexing barcode  (pos 3 to 4)
        rep_idx <- narrow(fq_R2read, 3, 4)
        barcode_string <-
            paste(sample_idx, umi_barcodes, rep_idx, sep = ":")
        # now trim the first 13 bp off the read and quality
        fq_R2read <- narrow(fq_R2read, 14, width(fq_R2))
        fq_R2qual <- narrow(fq_R2qual, 14, width(fq_R2))

    } else if (type == "RAMPAGE") {
        # trim barcodes first 15 bases from R2 and first 6 bases from R1 and prepare header
        # copy the index sequence (position 1 to 6 of R1)
        sample_idx <- narrow(fq_R1read, 1, 6)

        # copy pcr barcode (pos 1 to 15 of R2)
        umi_barcodes <- narrow(fq_R2read, 1, 15)
        barcode_string <- paste(sample_idx, umi_barcodes, sep = ":")

        # now trim the first 15 bp off the read and quality of R2
        fq_R2read <- narrow(fq_R2read, 16, width(fq_R2))
        fq_R2qual <- narrow(fq_R2qual, 16, width(fq_R2))
        fq_R1read <- narrow(fq_R1read, 7, width(fq_R1read))
        fq_R1qual <- narrow(fq_R1qual, 7, width(fq_R1qual))
    } else {
        message("type of protocol not MAPCap or RAMPAGE. Reads are not being trimmed.")
        barcode_string <- NA
    }

    outlist <- list(
        fq_R1read = fq_R1read,
        fq_R1qual = fq_R1qual,
        fq_R2read = fq_R2read,
        fq_R2qual = fq_R2qual,
        barcode_string = barcode_string,
        sample_idx = sample_idx
    )
    return(outlist)

}



#' Split paired-end fastq by barcodes
#'
#' @param expType character. experiment type (RAMPAGE, MAPCap or CAGE)
#' @param idx_name character. barcode ID
#' @param outfile_R1 character. output fastq file : Read 1
#' @param outfile_R2 character. output fastq file : Read 2
#' @param fastq_R1 character. input fastq file : Read 1
#' @param fastq_R2 character. input fastq file : Read 2
#' @param max_mismatch integer. max allowed mismatches
#'
#' @return kept reads corresponding to each barcode.
#'

split_fastq <- function(expType,
                        idx_name,
                        outfile_R1,
                        outfile_R2,
                        fastq_R1,
                        fastq_R2,
                        max_mismatch) {
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
        ## modify R1/R2 as per the protocol
        outlist <- get_newfastq(expType, fq_R1, fq_R2)

        ## make a new ShortReadQ object with new header and trimmed reads
        # new fastq R2
        fqid_R2 <- ShortRead::id(fq_R2)
        fqid_R2 <-
            sub(" ", "_", fqid_R2) # replace space by underscore in readName
        fq_R2new <- ShortRead::ShortReadQ(outlist$fq_R2read,
                                          outlist$fq_R2qual,
                                          Biostrings::BStringSet(paste(
                                                fqid_R2,
                                                outlist$barcode_string,
                                                sep = "#"
                                            )))

        # new fastq R1 (seq/qual not modified, just barcodes copied from R2)
        fqid_R1 <- ShortRead::id(fq_R1)
        fqid_R1 <-
            sub(" ", "_", fqid_R1) # replace space by underscore in readName

        fq_R1new <- ShortRead::ShortReadQ(outlist$fq_R1read,
                                          outlist$fq_R1qual,
                                          Biostrings::BStringSet(paste(
                                                fqid_R1,
                                                outlist$barcode_string,
                                                sep = "#"
                                            )))

        # demultiplex using given sample barcodes
        idx_name <- Biostrings::DNAString(idx_name)
        sample_idx <- Biostrings::DNAStringSet(outlist$sample_idx)
        id2keep <-
            as.logical(
                Biostrings::vcountPattern(idx_name,
                                            outlist$sample_idx,
                                            max.mismatch = max_mismatch)
            )

        #id2keep <- filter_byIDx(idx_name,
        #           fq_id = ShortRead::id(fq_R2),
        #           maxM = max_mismatch)
        # append to destination
        ShortRead::writeFastq(fq_R1new[id2keep], outfile_R1, "a")
        ShortRead::writeFastq(fq_R2new[id2keep], outfile_R2, "a")
        # add to count
        kept_reads <- kept_reads + sum(id2keep)
    }

    return(kept_reads)
}


#' Demultiplex and tag fastq files using sample barcodes
#' @rdname demultiplexFASTQ
#' @param CSobject CapSet object created using \code{\link{newCapSet}} function
#' @param outdir character. path to output directory
#' @param max_mismatch integer. maximum allowed mismatches in the sample barcode
#' @param ncores integrer. No. of cores/threads to use
#'
#' @return de-multiplxed fastq files corresponding to each barcode. The files are written
#'         on disk with the corresponding sample names as specified in the CapSet object
#'
#' @export
#' @examples
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#'
#' # demultiplex allowing one mismatch in sample indexes
#' dir.create("demult_fastq")
#' cs <- demultiplexFASTQ(cs, outdir =  "demult_fastq", max_mismatch = 1)
#'

setMethod("demultiplexFASTQ",
          signature = "CapSet",
          function(
                CSobject,
                outdir,
                max_mismatch,
                ncores) {
        protocol <- CSobject@expMethod
        sampleinfo <- sampleInfo(CSobject)
        destinations <- as.character(sampleinfo[, 1])
        idx_list <- as.character(rownames(sampleinfo))

        ## get the fastq to split (raise error if fastq untrimmed/not existing)
        fastq_R1 <- CSobject@fastq_R1
        fastq_R2 <- CSobject@fastq_R2

        # register parallel backend
        param <- getMCparams(ncores)
        if (!BiocParallel::bpisup(param)) {
        BiocParallel::bpstart(param)
        on.exit(BiocParallel::bpstop(param))
        }
        message("de-multiplexing the FASTQ file")
        ## filter and write
        info <-
            BiocParallel::bplapply(seq_along(destinations), function(i) {
                split1 <- file.path(outdir, paste0(destinations[i], "_R1.fastq.gz"))
                split2 <-
                    file.path(outdir, paste0(destinations[i], "_R2.fastq.gz"))
                ## stop if the outfiles already exist (otherwise the output would be appended)
                if (file.exists(split1) | file.exists(split2)) {
                    warning("Output files already exist! Overwriting..")
                    if (file.exists(split1))
                        file.remove(split1)
                    if (file.exists(split1))
                        file.remove(split2)
                }

                ## split files and save kept read number
                kept <- split_fastq(
                    expType = protocol,
                    idx_name = idx_list[i],
                    outfile_R1 = split1,
                    outfile_R2 = split2,
                    fastq_R1,
                    fastq_R2,
                    max_mismatch
                )
                dfout <- data.frame(R1 = split1,
                                    R2 = split2,
                                    kept_reads = kept)
                return(dfout)
            }, BPPARAM = param)

        ## add post-demult info to sampleInfo
        keptinfo <- do.call(rbind, info)
        sampleinfo$demult_reads <- keptinfo$kept_reads
        sampleinfo$demult_R1 <- as.character(keptinfo$R1)
        sampleinfo$demult_R2 <- as.character(keptinfo$R2)

        ## return object
        sampleInfo(CSobject) <- sampleinfo
        return(CSobject)

          }
)
