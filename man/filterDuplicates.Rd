% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/filterDuplicates.R
\name{filterDuplicates}
\alias{filterDuplicates}
\alias{filterDuplicates,CapSet-method}
\title{Filter PCR-duplicates from mapped files using internal UMIs}
\usage{
filterDuplicates(CSobject, outdir, ncores = 1, keepPairs = FALSE)

\S4method{filterDuplicates}{CapSet}(CSobject, outdir, ncores = 1, keepPairs = FALSE)
}
\arguments{
\item{CSobject}{an object of class \code{\link{CapSet}}}

\item{outdir}{character. output directory for filtered BAM files}

\item{ncores}{integer. No. of cores to use}

\item{keepPairs}{logical. indicating whether to keep pairs in the paired-end data.
(note: the pairs are treated as independent reads during duplicate removal).
Also use keepPairs = TRUE for single-end data.}
}
\value{
modified CapSet object with filtering information. Filtered BAM files are saved in `outdir`.
}
\description{
This script considers the read mapping start position and the UMI to determine whether a
             read is a PCR duplicate. All PCR duplicates are then removed and one entry per read is kept.
             In case of paired-end reads (MAPCap/RAMPAGE), only one end (R1) is kept after filtering, unless
             `keepPairs`` is set to TRUE
}
\examples{

# before running this
# 1. Create a CapSet object
# 2. de-multiplex the fastqs
# 3. map them

# load a previously saved CapSet object
cs <- exampleCSobject()
# filter duplicate reads from mapped BAM files
dir.create("filtered_bam")
cs <- filterDuplicates(cs, outdir = "filtered_bam")

}
