% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{readsTo3p}
\alias{readsTo3p}
\title{preprocess reads to count only 3' overlaps}
\usage{
readsTo3p(reads, width = 1, fix = "end")
}
\arguments{
\item{reads}{GAlignment object to resize}

\item{width}{integer. New read length}

\item{fix}{character. 'Start' for 5'}
}
\value{
Resized reads as GRanges
}
\description{
preprocess reads to count only 3' overlaps
}
