#' Calculate normalization factors from CapSet object
#'
#' @rdname getNormFactors
#' @param CSobject An object of class \code{\link{CapSet}}
#' @param features A \link[GenomicRanges]{GRanges-class}.object to count the reads on.
#' @param method Method to use for normalization. Options : "TMM","RLE","upperquartile","none"
#' @param ... Additional arguments passed to \link[edgeR]{calcNormFactors}
#'
#' @return Numeric vector of calculated normalization factors.
#' @importFrom GenomicAlignments summarizeOverlaps
#' @export
#'
#' @examples
#'
#'  # load a txdb object
#'  library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#'  seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
#'
#'  # get genes (only X chromsome, for simplicity)
#'  seqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "X"
#'  dm6genes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#'
#' # get norm factors by counting reads on genes
#' cs <- exampleCSobject()
#' normfacs <- getNormFactors(cs, dm6genes, method = "RLE")
#'

setMethod("getNormFactors",
        signature = "CapSet",
        function(CSobject,
                features,
                method,
                ...) {
        # get bam files
        si <- sampleInfo(CSobject)
        bam.files <- si$filtered_file

        ## get 5' read counts on the TSS from the bam.files
        counts <- summarizeOverlaps(features = features,
                                    reads = bam.files,
                                    preprocess.reads = ResizeReads)
        # make DGElist
        y <- edgeR::DGEList(counts = assay(counts))
        normfacs <- edgeR::calcNormFactors(y, method = method, ...)
        return(normfacs)
        }
)


#' @rdname fitDiffTSS
#' @param CSobject An object of class \code{\link{CapSet}}
#' @param TSSfile A .bed file with TSS positions to test for differential TSS analysis. If left
#'                empty, the union of detected TSS present within the provided CSobject would be plotted.
#' @param groups Character vector indicating the group into which each sample within the CSobject falls.
#'               the groups would be use to create a design matrix. As an example, replicates for one
#'               condition could be in the same group.
#' @param method Which method to use for differential expression analysis? options are "DESeq2" or "edgeR".
#'               If "DESeq2" is chosen, the library size is either estimated via DESeq2 (using "median of ratios")
#'               or can be provided via the "normFactors" option below. Setting the "normalization"
#'               (below) has no effect in that case.
#' @param normalization A character indicating the type of normalization to perform. Options are "windowTMM",
#'                 "TMM", "RLE", "upperquartile" or NULL (don't compute normalization factors).
#'                 If "windowTMM" is chosen, the normalization factors are calculated using the TMM
#'                 method on 10 kb windows of the genome. "TMM" computes TMM normalization using counts
#'                 from all the evaluated TSSs. If NULL, the external normalization factors can be used
#'                  (provided using `normFactors`).
#' @param normFactors external normalization factors (from Spike-Ins, for example).
#'
#' @param outplots Output pdf filename for plots. If provided, the plots for BCV, dispersion and
#'                 MDS plot is created and saved in this file.
#' @param plotRefSample Name of reference sample to plot for detection of composition bias in the
#'        data. Samples could be normalized using one of the provided normalization methods to
#'        control for composition bias.
#' @param ncores No. of cores/threads to use
#'
#' @return An object of class \link[edgeR]{DGEGLM-class} or
#'         \link[DESeq2]{DESeqDataSet-class}
#'
#' @importFrom graphics abline par plot smoothScatter
#' @importFrom grDevices dev.off pdf
#' @importFrom stats model.matrix
#' @importFrom methods as
#'
#' @export
#' @examples
#' # before running this
#' # 1. Create a CapSet object
#' # 2. de-multiplex the fastqs
#' # 3. map them
#' # 4. filter duplicate reads from mapped BAM
#' # 5. detect TSS
#' # 6. fit the diffTSS model
#' \dontrun{
#' # load a previously saved CapSet object
#' cs <- exampleCSobject()
#'
#' # count reads on all TSS (union) and fit a model using replicates within groups
#' csfit <- fitDiffTSS(cs, groups = rep(c("wt","mut"), each = 2), normalization = "internal",
#'                      outplots = NULL, plotRefSample = "embryo1")
#' save(csfit, file = "diffTSS_fit.Rdata")
#' }
#'

setMethod("fitDiffTSS",
          signature = "CapSet",
          function(CSobject,
             TSSfile,
             groups,
             method,
             normalization,
             normFactors,
             outplots,
             plotRefSample,
             ncores) {

        # get bam files and design
        si <- sampleInfo(CSobject)
        if (all(is.na(si$filtered_file))) {
            warning("Filtered files not found under sampleInfo(CSobject). Using mapped files")
            bam.files <- si$mapped_file
        } else {
            bam.files <- si$filtered_file
        }
        if (any(is.na(bam.files))) stop("Some or all of the bam files are not defined!")
        if (sum(file.exists(bam.files)) != length(bam.files)) {
            stop("One or more bam files don't exist! Check sampleInfo(CSobject) ")
        }

        samples <- as.character(si$samples)
        designDF <- data.frame(row.names = samples, group = groups)

        # Import tss locations to test
        if (is.null(TSSfile)) {
            if (is.null(CSobject@tss_detected)) stop("Detected TSS absent or not provided.")
            merged <- CSobject@tss_detected
            stopifnot(!(is.null(merged)))
            mergedall <- base::Reduce(S4Vectors::union, merged)
        } else {
            stopifnot(file.exists(TSSfile))
            mergedall <- rtracklayer::import.bed(TSSfile)
        }

        ## get 5' read counts on the TSS from the bam.files
        bp_param <- getMCparams(ncores)
        tsscounts <-
            GenomicAlignments::summarizeOverlaps(features = mergedall,
                                                 reads = bam.files,
                                                 singleEnd = TRUE,
                                                 preprocess.reads = ResizeReads,
                                                 BPPARAM = bp_param)
        ## check method
        if (method == "DESeq2") {
            counts <- assay(tsscounts)
            colnames(counts) <- samples
            dds <- DESeq2::DESeqDataSetFromMatrix(
                counts,
                rowRanges = SummarizedExperiment::rowRanges(tsscounts),
                colData = designDF,
                design = ~group)
            if (is.null(normFactors)) {
                message("Performing DESeq normalization")
                dds <- DESeq2::estimateSizeFactors(dds)
            } else {
                message("Using external normalization factors")
                DESeq2::sizeFactors(dds) <- normFactors
            }
            normalization <- "skip"
        } else if (method == "edgeR") {
            # make DGElist
            y <- edgeR::DGEList(counts = SummarizedExperiment::assay(tsscounts))
        }
        ## Get norm factors
        if (is.null(normalization)) {
            ## if normalization method not provided, look for external Norm factors
            normfacs <- normFactors
            # norm factors should be either NULL or a numeric vector
            if (!is.null(normfacs)) {
                stopifnot(is.numeric(normfacs) & length(normfacs) == length(groups))
            }
            y$samples$norm.factors <- normfacs

        } else if (normalization == "windowTMM") {
            ## normalization factors calculated using TMM on large Windows
            # fail early if plotRefSample is not given
            if (is.null(plotRefSample))
                stop("Please indicate reference sample for plotting of composition bias!")
            ## Internal normalization for composition bias : TMM
            # useful to try different bin sizes and see if the values are close to unity
            # (low composition effect)
            regionparam <- csaw::readParam(restrict = NULL,
                                           pe = ifelse(isTRUE(CSobject@paired_end), "first", "none"))
            binned <-
                csaw::windowCounts(bam.files,
                                   bin = TRUE,
                                   width = 10000,
                                   param = regionparam)
            normfacs <- csaw::normOffsets(binned) # close to unity
            names(normfacs) <- samples
            y.bin <- csaw::asDGEList(binned)

            y$samples$norm.factors <- normfacs
            if (!(is.na(plotRefSample))) {
            print(plotCompBias(dgelist = y.bin, samples, plotRefSample))
            }

        } else if (normalization %in% c("TMM", "RLE", "upperquartile", "none")) {
            ## normalization factors calculated using TMM on all TSSs
            y <- edgeR::calcNormFactors(y, method = normalization)
            if (!(is.na(plotRefSample))) {
                print(plotCompBias(dgelist = y, samples, plotRefSample))
            }

        } else if (normalization == "skip") {
            message("skipping additional normalizations")
        } else stop("Unknown normalization method!")

        if (method == "DESeq2") {
            ## proceed with DESeq2
            dds <- DESeq2::estimateDispersions(dds)
            fit <- DESeq2::nbinomWaldTest(dds)
            y <- NA

        } else if (method == "edgeR") {
            # proceed with edgeR
            designm <- model.matrix( ~ 0 + group, designDF)
            y <- edgeR::estimateDisp(y, designm)
            fit <- edgeR::glmQLFit(y, designm, robust = TRUE)
            # check prior degrees of freedom (avoid Infs)
            #message("Prior degrees of freedom : ")
            #print(summary(fit$df.prior))

        }

        ## make plots if asked
        if (!is.null(outplots)) {
            pdf(outplots)
            diffQCplots(method, fit, y, designMat = designDF)
            dev.off()
        }

        ## return the fit
        return(fit)
          }
    )


#' Make DESeq2 or edgeR QC plots
#' @importFrom ggplot2 ggplot aes aes_string geom_hline geom_vline geom_point labs theme_bw
#'                     scale_shape_manual
#' @param method one of "DESeq2" or "edgeR"
#' @param fit output of fitDiffTSS (if method = "edgeR")
#' @param y output of fitDiffTSS (if method = "DESeq2")
#'
#' @return Sparsity, dispersion and PCA plot (if method = DESeq2),
#'         BCV, dispersion and MDS plot (if method = "edgeR")
#'
diffQCplots <- function(method, fit, y, designMat) {

    if (method == "DESeq2") {

        ## plot sparsity
        DESeq2::plotSparsity(fit)
        ## plot Dispersions
        DESeq2::plotDispEsts(fit)
        ## plot PCA
        dds_rlog <- DESeq2::rlog(fit)
        PCAdata <- DESeq2::plotPCA(dds_rlog, intgroup = "group", returnData=TRUE)
        percentVar <- round(100 * attr(PCAdata, "percentVar"))
        print(ggplot(PCAdata, aes_string("PC1", "PC2", color = "group")) +
                  geom_hline(aes(yintercept = 0), colour = "grey") +
                  geom_vline(aes(xintercept = 0), colour = "grey") +
                  geom_point(size = 5) +
                  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
                       y = paste0("PC2: ", percentVar[2], "% variance"),
                       title = "PCA\n") +
                  theme_bw(base_size = 14) +
                  scale_shape_manual(values = c(0:18,33:17))
        )

    } else if (method == "edgeR") {

        ## plot BCV and fit
        par(mfrow = c(1, 2))
        o <- order(y$AveLogCPM)
        plot(
            y$AveLogCPM[o],
            sqrt(y$trended.dispersion[o]),
            type = "l",
            lwd = 2,
            ylim = c(0, 1),
            xlab = expression("Ave." ~ Log[2] ~ "CPM"),
            ylab = ("Biological coefficient of variation")
        )
        ## plot dispersion
        edgeR::plotQLDisp(fit)

        ## plot MDS (top genes/tss)
        par(mfrow = c(2, 2), mar = c(5, 4, 2, 2))
        for (top in c(100, 500, 1000, 5000)) {
            out <- limma::plotMDS(
                edgeR::cpm(y, log = TRUE),
                main = top,
                labels = designMat$group,
                top = top
            )
        }

    }
}


## for "windowTMM" norm : visualize Effect of TMM normalization on composition bias
plotCompBias <- function(dgelist, samples, plotRefSample) {

        normfacs <- dgelist$samples$norm.factors
        abundances <- edgeR::aveLogCPM(dgelist)
        adjc <- edgeR::cpm(dgelist, log = TRUE)
        colnames(adjc) <- samples

        # plot ref sample vs all other samples
        message("plotting the composition effect")
        sampnumber <- ncol(adjc) - 1
        cols_toplot <- grep(plotRefSample, samples, invert = TRUE)
        n <- ceiling(sampnumber / 3) #roundup to make divisible by 3

        par(cex.lab = 1.5, mfrow = c(n, 3))
        lapply(cols_toplot, function(x) {
            smoothScatter(
                abundances,
                adjc[, plotRefSample] - adjc[, x],
                ylim = c(-6, 6),
                xlab = "Average abundance",
                ylab = paste0("Log-ratio (", plotRefSample, " vs ", x, ")")
            )

        abline(h = log2(normfacs[plotRefSample] / normfacs[x]), col = "red")
        })
}

#' @rdname detectDiffTSS
#' @importFrom limma makeContrasts
#' @importFrom edgeR glmQLFTest
#' @importFrom ggplot2 ggplot aes_string geom_point geom_abline scale_color_manual labs theme_gray theme
#'
#' @export
#' @examples
#'
#' # before running this
#' # 1. Create a CapSet object
#' # 2. de-multiplex the fastqs
#' # 3. map them
#' # 4. filter duplicate reads from mapped BAM
#' # 5. detect TSS
#' # 6. fit the diff TSS model.
#'
#' \dontrun{
#' # load a previously saved DGEGLM object from step 5
#' csfit <- load("diffTSS_fit.Rdata")
#' dir <- system.file("extdata", package = "icetea")
#' # detect differentially expressed TSS between groups (return MA plot)
#' detectDiffTSS(csfit, testGroup = "mut", controlGroup = "wt",
#'                tssFile = file.path(dir, "testTSS_merged.bed"), MAplot_fdr = 0.05)
#'
#' }
#'

setMethod("detectDiffTSS",
          signature = "DGEGLM",
          function(fit,
            testGroup,
            contGroup,
            TSSfile,
            MAplot_fdr = NA) {
        # Import tss locations to test
        mergedall <- rtracklayer::import.bed(TSSfile)

        ## Testing the differential TSS
        # make contrast matrix
        contr <- paste0("group", testGroup , "-group", contGroup)
        contrast <-
            makeContrasts(contrasts = contr, levels = fit$design)
        # test
        results <- glmQLFTest(fit, contrast = contrast)
        top <- as.data.frame(edgeR::topTags(results, n = Inf))

        # sort output by pvalue and return
        difftss <- mergedall[as.numeric(rownames(top))]
        difftss$score <- top$FDR
        difftss$logFC <- top$logFC
        difftss$logCPM <- top$logCPM

        # MA plot
        if (!(is.na(MAplot_fdr))) {
            p <-
                ggplot(top, aes_string("logCPM", "logFC", col = factor(top$FDR < MAplot_fdr))) +
                geom_point(alpha = 0.5) +
                geom_abline(slope = 0, intercept = 0) +
                scale_color_manual(values = c("grey40", "darkred")) +
                labs(col = "Differentially Expressed") +
                theme_gray(base_size = 14) +
                theme(legend.position = "top")
            print(p)
        }

        return(difftss)
          }
)


#' @rdname detectDiffTSS
#' @export
#'
#' @importFrom DESeq2 results
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom ggplot2 ggplot aes_string geom_point geom_abline scale_color_manual labs
#'                     theme_gray theme scale_x_log10
#' @examples
#' \dontrun{
#' # load a previously saved DGEGLM object from step 5
#' csfit <- load("diffTSS_fit.Rdata")
#' dir <- system.file("extdata", package = "icetea")
#' # detect differentially expressed TSS between groups (return MA plot)
#' detectDiffTSS(csfit, testGroup = "mut", controlGroup = "wt", MAplot_fdr = 0.05)
#'
#' }
#'
setMethod("detectDiffTSS",
          signature = "DESeqDataSet",
          function(fit,
                   testGroup,
                   contGroup,
                   MAplot_fdr = NA) {

              # Get TSS location as GRanges from dds metadata
              mergedall <- rowRanges(fit)
              ## extract results from dds
              contrastvec <- c("group", testGroup, contGroup)
              ddr <- results(fit, contrast = contrastvec)

              # add information to TSS GRanges
              mergedall$score <- ddr$padj
              mergedall$logFC <- ddr$log2FoldChange
              mergedall$basemean <- ddr$baseMean

              # MA plot
              if (!(is.na(MAplot_fdr))) {
                  p <-
                      ggplot(as.data.frame(ddr), aes_string("baseMean", "log2FoldChange",
                                                            col = factor(ddr$padj < MAplot_fdr))) +
                      geom_point(alpha = 0.5) +
                      geom_abline(slope = 0, intercept = 0) +
                      scale_color_manual(values = c("grey40", "darkred")) +
                      labs(col = "Differentially Expressed") +
                      theme_gray(base_size = 14) +
                      theme(legend.position = "top") +
                      labs(x = "Mean Count", y = "log2(Fold-Change)") +
                      scale_x_log10()
                  print(p)
              }

              return(mergedall)
          }
)
