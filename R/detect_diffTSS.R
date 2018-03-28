#' Calculate normalization factors from CapSet object
#'
#' @rdname calcNormFactors
#' @param CSobject An object of class \code{\link{CapSet}}
#' @param features A \link[GenomicRanges]{GRanges}.object to count the reads on.
#' @param method Method to use for normalization. Options : "TMM","RLE","upperquartile","none"
#' @param ... Additional arguments passed to \link[edgeR]{calcNormFactors}
#'
#' @return Numeric vector of calculated normalization factors.
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
#' cs <- exampleCSObject()
#' normfacs <- calcNormFactors(cs, dm6genes, method = "RLE")
#'

setMethod("calcNormFactors",
        signature = "CapSet",
        function(CSobject,
                features,
                method,
                ...) {
        # get bam files
        si <- sampleInfo(CSobject)
        bam.files <- si$filtered_file

        ## get 5' read counts on the TSS from the bam.files
        counts <-
            GenomicAlignments::summarizeOverlaps(features = refRanges,
                                                 reads = bam.files,
                                                 preprocess.reads = ResizeReads)
        # make DGElist
        y <- edgeR::DGEList(counts = assay(tsscounts))
        normfacs <- edgeR::calcNormFactors(y, method = method, ...)
        return(normfacs)
        }
)


#' Detect differentially expressed Transcription Start Sites between two conditions (fit model)
#'
#' @rdname fitDiffTSS
#' @param CSobject An object of class \code{\link{CapSet}}
#' @param TSSfile A .bed file with TSS positions to test for differential TSS analysis. If left
#'                empty, the union of detected TSS present within the provided CSobject would be plotted.
#' @param groups Character vector indicating the group into which each sample within the CSobject falls.
#'               the groups would be use to create a design matrix. As an example, replicates for one
#'               condition could be in the same group.
#' @param normalization A character indicating the type of normalization to perform . Options are "windowTMM",
#'                 "globalTMM" or NULL (don't compute normalization factors).
#'                 If "windowTMM" is chosen, the normalization factors are calculated using the TMM
#'                 method on 10 kb windows of the genome. "globalTMM" compute TMM normalization using counts
#'                 from all the evaluated TSS. If NULL, the external normalization factors can be used (using
#'                 `normFactors`).
#' @param normFactors external normalization factors (from Spike-Ins, for example).
#'
#' @param outplots Output pdf filename for plots. If provided, the plots for BCV, dispersion and
#'                 MDS plot is created and saved in this file.
#' @param plotref Name of reference sample to plot for detection of composition bias in the
#'        data. Data is normalized using the TMM method to avoid composition bias.
#'
#' @return An object of class \link[edgeR]{DGEGLM}.
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
#'                      outplots = NULL, plotref = "embryo1")
#' save(csfit, file = "diffTSS_fit.Rdata")
#' }
#'

setMethod("fitDiffTSS",
        signature = "CapSet",
        function(CSobject,
            TSSfile,
            groups,
            normalization,
            normFactors,
            outplots,
            plotref) {

        # get bam files and design
        si <- sampleInfo(CSobject)
        bam.files <- si$filtered_file
        samples <- as.character(si$samples)
        design <- data.frame(row.names = samples, group = groups)

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
        tsscounts <-
            GenomicAlignments::summarizeOverlaps(features = mergedall,
                                                 reads = bam.files,
                                                 preprocess.reads = ResizeReads)
        # make DGElist
        y <- edgeR::DGEList(counts = assay(tsscounts))
        ## Get norm factors
        if (is.null(normalization)) {
            ## if normalization method not provided, look for external Norm factors
            normfacs <- normFactors
            # norm factors should be either NULL or a numeric vector
            if (!is.null(normfacs)) {
                stopifnot(is.numeric(normfacs) & length(normfacs) == length(groups))
            }
            y$samples$norm.factors <- normfacs
            plotCompBias(dgelist = y, samples, plotref)

        } else if (normalization == "windowTMM") {
            ## normalization factors calculated using TMM on large Windows
            # fail early if plotref is not given
            if (is.null(plotref))
                stop("Please indicate reference sample for plotting of composition bias!")
            ## Internal normalization for composition bias : TMM
            # useful to try different bin sizes and see if the values are close to unity (low composition effect)
            regionparam <- csaw::readParam(restrict = NULL)
            binned <-
                csaw::windowCounts(bam.files,
                                   bin = TRUE,
                                   width = 10000,
                                   param = regionparam)
            normfacs <- csaw::normOffsets(binned) # close to unity
            names(normfacs) <- samples
            y.bin <- csaw::asDGEList(binned)

            y$samples$norm.factors <- normfacs
            plotCompBias(dgelist = y.bin, samples, plotref)

        } else if (normalization %in% c("TMM", "RLE", "upperquartile", "none")) {
            ## normalization factors calculated using TMM on all TSSs
            y <- edgeR::calcNormFactors(y, method = normalization)
            plotCompBias(dgelist = y, samples, plotref)

        } else stop("Unknown normalization method!")

        # proceed with edgeR
        designm <- model.matrix( ~ 0 + group, design)
        y <- edgeR::estimateDisp(y, designm)
        fit <- edgeR::glmQLFit(y, designm, robust = TRUE)
        # check prior degrees of freedom (avoid Infs)
        message("Prior degrees of freedom : ")
        print(summary(fit$df.prior))

        ## make plots if asked
        if (!is.null(outplots)) {
            pdf(outplots)
            ## check that the Fit is good
            par(mfrow = c(1, 2))
            # BCV plot
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
            # dispersion plot
            edgeR::plotQLDisp(fit)

            ## check with MDSplot if there is batch effect
            par(mfrow = c(2, 2), mar = c(5, 4, 2, 2))
            for (top in c(100, 500, 1000, 5000)) {
                out <- limma::plotMDS(
                    edgeR::cpm(y, log = TRUE),
                    main = top,
                    labels = design$group,
                    top = top
                )
            }

            dev.off()
        }

        ## return the fit
        return(fit)

        }
)

## visualize Effect of TMM normalization on composition bias
plotCompBias <- function(dgelist, samples, plotref) {

        normfacs <- y$samples$norm.factors
        abundances <- edgeR::aveLogCPM(dgelist)
        adjc <- edgeR::cpm(dgelist, log = TRUE)
        colnames(adjc) <- samples

        # plot ref sample vs all other samples
        message("plotting the composition effect")
        sampnumber <- ncol(adjc) - 1
        cols_toplot <- grep(plotref, samples, invert = TRUE)
        n <- ceiling(sampnumber / 3) #roundup to make divisible by 3

        par(cex.lab = 1.5, mfrow = c(n, 3))
        lapply(cols_toplot, function(x) {
            smoothScatter(
                abundances,
                adjc[, plotref] - adjc[, x],
                ylim = c(-6, 6),
                xlab = "Average abundance",
                ylab = paste0("Log-ratio (", plotref, " vs ", x, ")")
            )

        abline(h = log2(normfacs[plotref] / normfacs[x]), col = "red")
        })
}

#' Detect differentially expressed Transcription Start Sites between two conditions (test)
#'
#' @rdname detectDiffTSS
#' @param fit DGEGLM object (output of \code{\link{fitDiffTSS}} command )
#' @param testGroup Test group name
#' @param contGroup Control group name
#' @param TSSfile The TSS .bed file used for \code{\link{fitDiffTSS}} command
#' @param MAplot_fdr FDR threshold to mark differentially expressed TSS in MAplot (NA = Don't make an MAplot)
#'
#' @return A \code{\link{GRanges}} object containing p-values of differential expression for each TSS.
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
