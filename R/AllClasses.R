#' Create a new CapSet object
#'
#' The function creates an object of class `CapSet`, used for the TSS analysis.
#' A CapSet object can be created using the the raw, multiplxed fastq files along
#' with the list of sample indexes and corresponding sample names. In case the
#' files are already de-multiplexed or mapped, the CapSet object can also be created
#' using the path to demultiplexed fastq/mapped or filtered BAM files, along with
#' corresponding sample names. In these cases statistics and operations for the
#' missing files would not be possible.
#'
#'
#' @param expMethod experiment method ('CAGE', 'RAMPAGE' or 'MAPCap')
#' @param fastq_R1 path for Read R1 (or file path for single end reads)
#' @param fastq_R2 path for Read R2 (for paired end reads)
#' @param idxList a vector of index sequences (for demultiplexing)
#' @param sampleNames (required) a vector of sample names corresponding to the provided files
#' @param demult_R1 a vector of file paths for demultiplexed R1 reads
#' @param demult_R2 a vector of file paths for demultiplexed R2 reads
#' @param mapped_file a vector of file paths for mapped BAM files.
#' @param filtered_file a vector of file paths for de-duplicated BAM files.
#' @param paired_end logical, indiciting whether the data is paired end
#'
#' @return An object of class CapSet
#' @importFrom methods new is
#' @importFrom Rsamtools countBam ScanBamParam scanBamFlag BamFileList
#' @export
#'
#' @examples
#'
#' # list of barcode IDs
#' idxlist <- c("CAAGTG", "TTAGCC", "GTGGAA", "TGTGAG")
#'
#' dir <- system.file("extdata", package="icetea")
#' # corresponding sample names
#' fnames <- c("embryo1", "embryo2", "embryo3", "embryo4")
#'
#' ## CapSet object from raw (multiplexed) fastq files
#' cs <- newCapSet(expMethod = 'MAPCap',
#'        fastq_R1 = file.path(dir, 'mapcap_test_R1.fastq.gz'),
#'        fastq_R2 = file.path(dir, 'mapcap_test_R2.fastq.gz'),
#'        idxList = idxlist,
#'        sampleNames = fnames)
#'
#' ## CapSet object from mapped BAM files
#' bams <- list.files(file.path(dir, 'bam'), pattern = '.bam$', full.names = TRUE)
#' cs <- newCapSet(expMethod = 'MAPCap',
#'        mapped_file = bams,
#'        sampleNames = fnames)
#'

newCapSet <- function(
                    expMethod,
                    fastq_R1 = NULL,
                    fastq_R2 = NULL,
                    idxList = NULL,
                    sampleNames,
                    demult_R1 = NA,
                    demult_R2 = NA,
                    mapped_file = NA,
                    filtered_file = NA,
                    paired_end = TRUE) {
    ## Get numbers from already provided files
    suppressWarnings({

        if (!is.null(idxList)) {
            # Check validity of index seq
            idxList <- as.character(Biostrings::DNAStringSet(idxList))
            # the above func would throw an error if seq not valid
        }

        # R1 (read only if none of the entries are NA)
        if (!any(is.na(demult_R1))) {
            message("Checking de-multiplexed R1 reads")
            if (sum(vapply(demult_R1, file.exists, logical(1L))) != length(demult_R1)) {
                stop("One or more R1 read files don't exist!")
            }
            r1_counts <- as.integer(ShortRead::countLines(demult_R1)/4)
        } else {
            r1_counts <- NA
        }
        # R2 (read only if none of the entries are NA)
        if (!any(is.na(demult_R2))) {
            message("Checking de-multiplexed R2 reads")
            if (sum(vapply(demult_R2, file.exists, logical(1L))) != length(demult_R2)) {
                stop("One or more R2 read files don't exist!")
            }

            r2_counts <- as.integer(ShortRead::countLines(demult_R2)/4)
            # R1 and R2 have same no of reads?
            if (!all(r2_counts == r1_counts)) {
                stop("discrepency between number of R2 and R1 reads in fastq files!")
            }

        } else {
            # paired_end stays FALSE
            r2_counts <- NA
        }
        # BAM (read only if none of the entries are NA)
        if (!any(is.na(mapped_file))) {
            message("Checking mapped file")
            if (sum(vapply(mapped_file, file.exists, logical(1L))) != length(mapped_file)) {
                stop("One or more mapped files don't exist!")
            }

            ## count all reads if single end, else count only R1
            mapped_readcounts <- countBam(BamFileList(mapped_file),
                                            param = ScanBamParam(
                                              flag = getBamFlags(countAll = !paired_end)
                                            ))[, 6] # "records"
        } else {
            mapped_readcounts <- NA
        }
        # Filtered BAM (read only if none of the entries are NA)
        if (!any(is.na(filtered_file))) {
            message("Checking de-duplicated file")
            if (sum(vapply(filtered_file, file.exists, logical(1L))) != length(filtered_file)) {
                stop("One or more de-duplicated files don't exist!")
            }
            ## count all reads if single end, else count only R1
            filt_readcounts <- countBam(BamFileList(filtered_file),
                                        param = ScanBamParam(
                                            flag = getBamFlags(countAll = !paired_end)
                                        ))[, 6] # "records"
        } else {
            filt_readcounts <- NA
        }

    })

    # make sampleinfo DataFrame
    info <- S4Vectors::DataFrame(
        row.names = idxList,
        samples = sampleNames,
        demult_R1 = demult_R1,
        demult_R2 = demult_R2,
        mapped_file = mapped_file,
        filtered_file = filtered_file,
        demult_reads = r1_counts,
        num_mapped = mapped_readcounts,
        num_filtered = filt_readcounts,
        num_intss = NA
    )

    # get fastq type (single or paired)
    #fastqType <- ifelse(is.null(fastq_R2), "single", "paired")
    # create an instance of CapSet
    new(
        "CapSet",
        sampleInfo = info,
        # sample Information
        paired_end = paired_end,
        fastq_R1 = fastq_R1,
        fastq_R2 = fastq_R2,
        expMethod = expMethod,
        tss_detected = NULL
    )
}


#' Check capset validity
#'
#' @param object capset object
#'
#' @return boolean
#'

check_capSet <- function(object) {
    errors <- character()

    ## extract slots
    fqtype <- object@fastqType
    R1 <- object@fastq_R1
    R2 <- object@fastq_R2
    exp <- object@expMethod
    info <- object@sampleInfo
    tss <- object@tss_detected

    ## validate slots
    # experiment
    if (!(exp %in% c("CAGE", "RAMPAGE", "MAPCap"))) {
        msg <-
            paste0("Experiment type should be among : 'CAGE', 'RAMPAGE' or 'MAPCap' ")
        errors <- c(errors, msg)
    }
    # sampleInfo
    if (!is(info, "DataFrame")) {
        msg <- paste0("sampleInfo should be a DataFrame object ")
        errors <- c(errors, msg)
    }
    si_names <- c("samples", "demult_R1", "demult_R2",
                    "mapped_file", "filtered_file", "demult_reads",
                    "num_mapped", "num_filtered", "num_intss" )
    if (!(all(colnames(info) %in% si_names))) {
        stop("Column names in sampleInfo DataFrame are not correct!")
    }
    si_types <- c("character", "integer", "numeric", "logical")
    col_classes <- vapply(info, class, character(1L))
    if (!(all(col_classes %in% si_types))) {
        stop("Column classes in sampleInfo DataFrame are not correct!")
    }
    # TSS info
    if (!(class(tss) %in% c("NULL", "CompressedGRangesList", "GRangesList"))) {
        msg <- paste0("tss_detected should be a GRangesList object ")
        errors <- c(errors, msg)
    }

    ## return
    if (length(errors) == 0)
        TRUE
    else
        errors
}

## char or NULL class
setClassUnion("charOrNULL", c("character", "NULL"))

#' CapSet object
#'
#' @rdname newCapSet
#' @slot fastqType Type of fastq ('single' or 'paired')
#' @slot fastq_R1 Path to R1 fastq
#' @slot fastq_R2 Path to R1 fastq (for paired-end data)
#' @slot expMethod Name of protocol (RAMPGE or MAPCap)
#' @slot sampleInfo A DataFrame object created using information from
#'                  \code{\link{newCapSet}} function
#' @slot tss_detected A GRangesList object of detected TSS
#' @importClassesFrom S4Vectors DataFrame
#'
CapSet <- setClass(
    "CapSet",
    slots = c(
        fastqType = "character",
        fastq_R1 = "charOrNULL",
        fastq_R2 = "charOrNULL",
        expMethod = "character",
        sampleInfo = "DataFrame",
        paired_end = "logical",
        tss_detected = "ANY"
    ),
    validity = check_capSet
)
