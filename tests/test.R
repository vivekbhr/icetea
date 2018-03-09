#dir <- system.file("extdata", package = "icetea")
#si <- data.frame(row.names = c("CAAGTG", "TTAGCC", "GTGGAA", "TGTGAG"),
#		 fnames = c("embryo1", "embryo2", "embryo3", "embryo4"))
## working examples
#cs <- newCapSet(expMethod = 'MAPCap', fastqType = 'paired',
#		fastq_R1 = file.path(dir, 'mapcap_test_R1.fastq.gz'),
#		fastq_R2 = file.path(dir, 'mapcap_test_R2.fastq.gz'),
#		sampleInfo = si)
#cs <- demultiplex_fastq(cs, max_mismatch = 1, outdir = file.path(dir, 'demult_fastq'))
#cs <- mapCaps(cs, "/data/akhtar/bhardwaj/my_annotations/drosophila_dm6/subread_index/dm6",
#	      file.path(dir, "bam"), nthreads = 10)
#cs <- filterDuplicates(cs, file.path(dir, "filtered_bam"))
#cs <- detect_TSS(cs, groups = c("wt", "wt", "mut", "mut"),
#		 restrictChr = 'X', foldChange = 4)
#export_tss(cs, outfile_prefix = "testTSS")

#library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
#library(GenomicFeatures)
#seqlevelsStyle(dm6trans_all) <- "ENSEMBL"

#dm6trans_bygene <- transcriptsBy(TxDb.Dmelanogaster.UCSC.dm6.ensGene, "gene")
#dm6trans_all <- transcripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

#gcounts <- get_geneCounts(cs, dm6trans)
#csfit <- fitDiffTSS(cs, "testTSS_merged.bed", groups = c("wt", "wt", "mut", "mut"),
#		     outplots = "test", plotref = "embryo1")

# plots
#plot_TSSprecision(reference = dm6trans_all, detectedTSS = cs)
#plot_readStats(cs)

