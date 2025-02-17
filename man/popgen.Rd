% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popgen.R
\name{popgen}
\alias{popgen}
\title{Population Generation}
\usage{
popgen(n, p)
}
\arguments{
\item{n}{number of genes in the chromosome associated with a solution}

\item{p}{number of elements (chromosomes) in the population}
}
\value{
\item{pop}{ = p x n real-valued matrix}
}
\description{
Generates p vectors of random-keys
}
\details{
popgen
}
\examples{
X<-popgen(n=20,p=10)
}
\references{
1. Gonçalves, J.R. and Resende, M.G.C. (2011). Biased random-key genetic algorithms for combinatorial optimization,
}
