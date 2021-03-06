#' Map the data from 5' profiling techniques
#'
#' @rdname mapCaps
#' @param CSobject An object of class \code{\link{CapSet}}
#' @param genomeIndex character. Path to the Subread index file. Should end with the basename of the index.
#' @param outdir character. Output directory path
#' @param ncores integer. Number of cores/threads to use for mapping.
#' @param logfile character. A log file to write the processing message.
#' @param externalGTF character. provide external annotation file in `GTF` format , if present
#'                    to increase alignment accuracy
#'
#' @return modified CapSet object with mapping information. Mapped and sorted BAM
#'         files are saved in `outdir`.
#'
#' @importFrom utils capture.output
#' @importFrom methods validObject is
#' @importFrom Rsamtools sortBam indexBam countBam ScanBamParam scanBamFlag
#'
#' @examples
#' \dontrun{
#' # before mapping :
#' # 1. Create a CapSet object
#' # 2. de-multiplex the fastqs
#'
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#'
#' # map the data (not available on windows)
#' library(Rsubread)
#' dir.create("bam")
#' buildindex(basename = "dm6", reference = "/path/to/dm6_genome.fa")
#' cs <- mapCaps(cs, genomeIndex = "dm6", outdir = "bam", nthreads = 10)
#'
#' }
#'

setMethod("mapCaps",
        signature = "CapSet",
            function(
                CSobject,
                genomeIndex,
                outdir,
                externalGTF,
                ncores,
                logfile) {
        ## extract info
        sampleInfo <- sampleInfo(CSobject)
        expMethod <- CSobject@expMethod

        # test whether the CapSet object has demultiplexed fastqs
        demult <- !(is.null(sampleInfo$demult_R1)) # TRUE/FALSE

        if (demult) {
            if (sum(vapply(sampleInfo$demult_R1, file.exists, logical(1))) != nrow(sampleInfo)) {
                stop(
                    "One or more demultiplxed fastq files don't exist.
                    Please check file paths under sampleInfo(CSobject)"
                )
            }
            }

        if (demult) {
            samplelist <- sampleInfo$samples
            R1_list <- as.character(sampleInfo$demult_R1)
            R2_list <- as.character(sampleInfo$demult_R2)
        } else {
            samplelist <- list("raw")
            R1_list <- CSobject@fastq_R1
            R2_list <- CSobject@fastq_R2
        }
        # get annotation if provided
        if (!is.null(externalGTF)) {
            useAnnot <- TRUE
        } else {
            useAnnot <- FALSE
        }
        mapstat <- mapply(function(sample, R1, R2) {
            message(paste0("Mapping sample : ", sample))
            # If R2 = NA, convert it to NULL such that subread uses single-end
            if (is.na(R2)) R2 <- NULL
            # Align using RSubread
            tmpout <- file.path(outdir, paste0(sample, ".tmp.bam"))
            capture.output(
                Rsubread::subjunc(
                    index = genomeIndex,
                    # change here
                    readfile1 = R1,
                    readfile2 = R2,
                    output_file = tmpout,
                    annot.ext = externalGTF,
                    isGTF = TRUE,
                    # annotation (v1.28 onwards)
                    useAnnotation = useAnnot,
                    annot.inbuilt = NULL,
                    # neglected
                    GTF.featureType = "exon",
                    GTF.attrType = "gene_id",
                    chrAliases = NULL,
                    input_format = "gzFASTQ",
                    output_format = "BAM",
                    nsubreads = 20,
                    TH1 = 1,
                    TH2 = 1,
                    maxMismatches = 3,
                    nthreads = ncores,
                    indels = 5,
                    complexIndels = FALSE,
                    phredOffset = 33,
                    unique = TRUE,
                    nBestLocations = 1,
                    minFragLength = 10,
                    maxFragLength = 600,
                    PE_orientation = "fr",
                    nTrim5 = 0,
                    nTrim3 = 0,
                    readGroupID = NULL,
                    readGroup = NULL,
                    color2base = FALSE,
                    # subjunc-specific
                    reportAllJunctions = FALSE,
                ),
                file = logfile,
                append = TRUE
            )

            # Sort and Index
            message("Sorting and Indexing")
            dest <- file.path(outdir, sample)
            sortBam(file = tmpout, destination = dest) # adds .bam suffix
            indexBam(paste0(dest, ".bam"))
            file.remove(tmpout)

            # Get mapping stats
            stat <- countBam(paste0(dest, ".bam"),
                            param = ScanBamParam(
                                flag = scanBamFlag(
                                    isUnmappedQuery = FALSE,
                                    isFirstMateRead = TRUE,
                                    isSecondaryAlignment = FALSE
                                 )
                            ))[, 5:6] # "file" and "records"
            stat$file <- as.character(stat$file)
            stat$records <- as.integer(stat$records)
            return(stat)
        }, samplelist, R1_list, R2_list)

        # edit sampleinfo of CapSet
        si <- sampleInfo(CSobject)
        maptable <- as.data.frame(t(mapstat))
        si$mapped_file <- file.path(outdir, as.character(maptable$file))
        si$num_mapped <- as.integer(maptable$records)
        sampleInfo(CSobject) <- si

        validObject(CSobject)
        return(CSobject)
        }
)
