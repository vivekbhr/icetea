context("icetea")

test_that("De-duplication results in expected No. of reads", {
    cs <- exampleCSobject()
    cs2 <- filterDuplicates(cs, outdir = ".", keepPairs = FALSE)
    # producs non-empty files?
    filtbams <- sampleInfo(cs2)$filtered_file
    expect_true(all(file.exists(filtbams)))
    # expected No. of reads?
    expect_equal(sampleInfo(cs)$num_filtered, sampleInfo(cs2)$num_filtered)
})
