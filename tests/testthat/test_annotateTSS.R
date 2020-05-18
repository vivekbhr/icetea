context("icetea")

test_that("TSS annotation produces correct table", {
    library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
    seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
    seqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "X"
    dir <- system.file("extdata", package = "icetea")
    tss <- file.path(dir, "testTSS_merged.bed")
    annot <- annotateTSS(tss, txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
                         plotValue = "number", outFile = "TSS_annot.pdf")
    # expect right output
    expect_is(annot, "data.frame")
    expected_vals <- c(0, 7, 56, 3, 4, 3, 8)
    expect_equal(annot$value, expected_vals)

    })
