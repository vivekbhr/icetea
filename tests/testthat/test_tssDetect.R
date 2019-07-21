context("icetea")

test_that("TSS calling produces expected No. of reads and peaks", {
    cs <- exampleCSobject()
    cs2 <- detectTSS(cs, groups = rep(c("wt","mut"), each = 2), outfile_prefix = "testTSS",
            foldChange = 2, restrictChr = "X", ncores = 1, sliding = FALSE)
    # producs .Rdata file?
    expect_true(file.exists("testTSS.Rdata"))

    # expected No. of peaks in each group?
    expect_equal(length(cs2@tss_detected), 2)
    vec <- c(1,2)
    names(vec) <- c("wt", "mut")
    expect_identical(vapply(cs2@tss_detected, length, numeric(1L)), vec)

    # expected No. of reads in peaks?
    expect_equal(sampleInfo(cs)$num_intss, sampleInfo(cs2)$num_intss)
})
