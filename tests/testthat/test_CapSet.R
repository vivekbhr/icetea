context("icetea")

test_that("CapSet object produces correct errors", {

    ## correct entries
    dir = system.file("extdata", package = "icetea")
    index = c("CAAGTG", "TTAGCC", "GTGGAA", "TGTGAG")
    fnames = c("embryo1", "embryo2", "embryo3", "embryo4")
    exp = "MAPCap"
    r1 = file.path(dir, 'mapcap_test_R1.fastq.gz')
    r2 = file.path(dir, 'mapcap_test_R2.fastq.gz')
    ## passing example
    cs <- newCapSet(exp, r1, r2,index, fnames)
    expect_type(cs, "S4")

    ## failing examples
    exp2 = "test"
    expect_error(newCapSet(exp2, r1, r2, index, fnames))
    fnames2 <- fnames[-1]
    expect_error(newCapSet(exp, r1, r2, index, fnames2))
    index2 <- index
    index2[1] <- "ATCGXX"
    expect_error(newCapSet(exp, r1, r2, index2, fnames))
    ## sampleInfo setter and getter
    si <- sampleInfo(cs)
    expect_type(si, "DataFrame")
    si$num_intss <- as.matrix(si$num_intss)
    expect_error(sampleInfo(cs) <- si)
})

