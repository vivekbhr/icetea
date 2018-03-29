context("icetea")

test_that("De-multiplexing results in expected No. of reads", {
    cs <- exampleCSobject()
    cs2 <- demultiplexFASTQ(cs, outdir = ".", max_mismatch = 1)
    # producs non-empty files?
    r1 <- sampleInfo(cs2)$demult_R1
    r2 <- sampleInfo(cs2)$demult_R2
    expect_true(all(file.exists(r1)))
    expect_true(all(file.exists(r2)))
    # expected No. of reads?
    expect_equal(sampleInfo(cs)$demult_reads, sampleInfo(cs2)$demult_reads)
})
