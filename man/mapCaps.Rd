% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/mapCaps.R
\name{mapCaps}
\alias{mapCaps}
\alias{mapCaps,CapSet-method}
\title{Map the data from 5' profiling techniques}
\usage{
mapCaps(
  CSobject,
  genomeIndex,
  outdir,
  externalGTF = NULL,
  ncores = 1,
  logfile = NULL
)

\S4method{mapCaps}{CapSet}(
  CSobject,
  genomeIndex,
  outdir,
  externalGTF = NULL,
  ncores = 1,
  logfile = NULL
)
}
\arguments{
\item{CSobject}{An object of class \code{\link{CapSet}}}

\item{genomeIndex}{character. Path to the Subread index file. Should end with the basename of the index.}

\item{outdir}{character. Output directory path}

\item{externalGTF}{character. provide external annotation file in `GTF` format , if present
to increase alignment accuracy}

\item{ncores}{integer. Number of cores/threads to use for mapping.}

\item{logfile}{character. A log file to write the processing message.}
}
\value{
modified CapSet object with mapping information. Mapped and sorted BAM
        files are saved in `outdir`.
}
\description{
Map the data from 5' profiling techniques
}
\examples{
\dontrun{
# before mapping :
# 1. Create a CapSet object
# 2. de-multiplex the fastqs

# load a previously saved CapSet object
cs <- exampleCSobject()

# map the data (not available on windows)
library(Rsubread)
dir.create("bam")
buildindex(basename = "dm6", reference = "/path/to/dm6_genome.fa")
cs <- mapCaps(cs, genomeIndex = "dm6", outdir = "bam", nthreads = 10)

}

}
