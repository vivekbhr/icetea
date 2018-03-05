#' Detect differentially expressed Transcription Start Sites between two conditions (fit model)
#'
#' @param CSobject An object of class \code{\link{CapSet}}
#' @param TSSfile A .bed file with TSS positions to test for differential TSS analysis. If left
#'                empty, the union of detected TSS present within the provided CSobject would be plotted.
#' @param groups Character vector indicating the group into which each sample within the CSobject falls.
#'               the groups would be use to create a design matrix. As an example, replicates for one
#'               condition could be in the same group.
#' @param normalization Either a character indicating the type of normalization to perform ("internal",
#'                 "external" or "none"), or a numeric vector vector with pre-computed normalization factors.
#'                 If internal normalization is chosen, the normalization factors are calculated using the TMM
#'                 method on large windows of the genome. For "external" normalization, the normalization factors
#'                 from the provided spike-in samples are used. "none" performs no normalization.
#'
#' @param CSobjectSpikeIn Another CapSet object produced using the spike-in mapping.
#'
#' @param outplots Output pdf filename for plots. If provided, the plots for BCV, dispersion and
#'                 MDS plot is created and saved in this file.
#' @param plotref Name of reference sample to plot for detection of composition bias in the
#' 		  data. Data is normalized using the TMM method to avoid composition bias.
#'
#' @return Returns an object of class DGEGLM.
#'
#' @export
#' @importFrom graphics abline par plot smoothScatter
#' @importFrom grDevices dev.off pdf
#' @importFrom stats model.matrix
#' @importFrom methods as
#'
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
#' csfit <- fit_diffTSS(cs, groups = rep(c("wt","mut"), each = 2), normalization = "internal",
#'                      outplots = NULL, plotref = "embryo1")
#' save(csfit, file = "diffTSS_fit.Rdata")
#' }
#'

fit_diffTSS <-
    function(CSobject,
             TSSfile = NULL,
             groups,
             normalization = "internal",
             CSobjectSpikeIn,
             outplots = NULL,
             plotref) {
        ## assert the input
        stopifnot(is(CSobject, "CapSet"))

        # get bam files and design
        si <- sampleInfo(CSobject)
        bam.files <- si$filtered_file
        samples <- as.character(si$samples)
        design <- data.frame(row.names = samples, group = groups)

        # Import tss locations to test
        if (is.null(TSSfile)) {
            merged <- CSobject@tss_detected
            stopifnot(!(is.null(merged)))
            mergedall <- base::Reduce(S4Vectors::union, merged)
        } else {
            mergedall <- rtracklayer::import.bed(TSSfile)
        }

        if (normalization == "internal") {
            # fail early if plotref is not given
            if (is.null(plotref))
                stop("Please indicate reference sample for plotting of composition bias!")
            ## Internal normalization for composition bias : TMM
            # useful to try different bin sizes and see if the values are close to unity (low composition effect)
            regionparam <- csaw::readParam(minq = 30, restrict = NULL)
            binned <-
                csaw::windowCounts(bam.files,
                                   bin = TRUE,
                                   width = 10000,
                                   param = regionparam)
            normfacs <- csaw::normOffsets(binned) # close to unity
            names(normfacs) <- samples

            ## visualize Effect of TMM normalization on composition bias
            y.bin <- csaw::asDGEList(binned)
            bin.ab <- edgeR::aveLogCPM(y.bin)
            adjc <- edgeR::cpm(y.bin, log = TRUE)
            colnames(adjc) <- samples

            # plot ref sample vs all other samples
            message("plotting the composition effect")
            sampnumber <- ncol(adjc) - 1
            cols_toplot <- grep(plotref, samples, invert = TRUE)
            n <- ceiling(sampnumber / 3) #roundup to make divisible by 3

            par(cex.lab = 1.5, mfrow = c(n, 3))
            lapply(cols_toplot, function(x) {
                smoothScatter(
                    bin.ab,
                    adjc[, plotref] - adjc[, x],
                    ylim = c(-6, 6),
                    xlab = "Average abundance",
                    ylab = paste0("Log-ratio (", plotref, " vs ", x, ")")
                )

                abline(h = log2(normfacs[plotref] / normfacs[x]),
                       col = "red")

            })

        } else {
            normfacs <- NULL
        }

        ## get 5' read counts on the locations from the bam.files
        # function to resize reads
        ResizeReads <- function(reads,
                                width = 1,
                                fix = "start",
                                ...) {
            reads <- as(reads, "GRanges")
            stopifnot(all(GenomicRanges::strand(reads) != "*"))
            GenomicRanges::resize(reads, width = width, fix = fix, ...)
        }

        # now read the data
        tsscounts <-
            GenomicAlignments::summarizeOverlaps(features = mergedall,
                                                 reads = bam.files,
                                                 preprocess.reads = ResizeReads)
        #### ------ Now do EdgeR ------ ####
        y <- csaw::asDGEList(tsscounts, norm.factors = normfacs)
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



#' Detect differentially expressed Transcription Start Sites between two conditions (test)
#'
#' @param fit DGEGLM object (output of \code{\link{fit_diffTSS}} command )
#' @param testGroup Test group name
#' @param contGroup Control group name
#' @param TSSfile The TSS .bed file used for \code{\link{fit_diffTSS}} command
#' @param MAplot_fdr FDR threshold to mark differentially expressed TSS in MAplot (NA = Don't make an MAplot)
#'
#' @return A \code{\link{GRanges}} object containing p-values of differential expression for each TSS.
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_point geom_abline scale_color_manual labs theme_gray theme
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
#' detect_diffTSS(csfit, testGroup = "mut", controlGroup = "wt",
#'                tssFile = file.path(dir, "testTSS_merged.bed"), MAplot_fdr = 0.05)
#'
#' }
#'

detect_diffTSS <-
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
            limma::makeContrasts(contrasts = contr, levels = fit$design)
        # test
        results <- edgeR::glmQLFTest(fit, contrast = contrast)
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
