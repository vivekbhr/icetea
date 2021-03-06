% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateTSS.R
\name{annotateTSS}
\alias{annotateTSS}
\title{Annotate the provided Transcription Start Sites}
\usage{
annotateTSS(
  tssBED,
  txdb,
  featureRank = c("fiveUTR", "promoter", "intron", "coding", "spliceSite", "threeUTR",
    "intergenic"),
  plotValue = "number",
  outFile = NULL
)
}
\arguments{
\item{tssBED}{A bed file with detected TSS/differential TSS coordinates}

\item{txdb}{A txdb object.}

\item{featureRank}{A vector with features to use for breaking ties, in decending order of
preference (highest to lowest),}

\item{plotValue}{What values to plot (choose from "number", "percent" or NULL for no plot)}

\item{outFile}{Output file name. (filename extention would be used to determine type).
If outfile not specified, the plot would be retured on the screen}
}
\value{
A data.frame with number of TSS falling into each feature
}
\description{
This function annotates the provided TSS bed file to provide the number of TSS
falling within the genomic features from a given TxDB object. In order to break
ties between overlapping features, the function ranks the features by preference.
By default, the following order is used: fiveUTR > promoter > intron > coding >
spliceSite > threeUTR > intergenic. A custom order of feature ranks can also be
provided.
}
\examples{
# load a txdb object
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "ENSEMBL"
# limiting the annotation to X chromosome
seqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "X"

# annotate a given TSS bed file
dir <- system.file("extdata", package = "icetea")
tss <- file.path(dir, "testTSS_merged.bed")
annotations <- annotateTSS(tssBED = tss, TxDb.Dmelanogaster.UCSC.dm6.ensGene,
               plotValue = "number", outFile = "TSS_annot.pdf")

}
