% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mungeCARdata4stan_irregular.R
\name{mungeCARdata4stan_irregular}
\alias{mungeCARdata4stan_irregular}
\title{Prepare necessary components of the IAR model for use in \code{stan}.}
\usage{
mungeCARdata4stan_irregular(adjBUGS, numBUGS)
}
\arguments{
\item{adjBUGS}{A vector of neighbor indices.}

\item{numBUGS}{A vector containing the number of neighbors for each point.}
}
\value{
Returns a list where \code{J} is the number predictors,
\code{J_edges} is the number of edges in connected graph representing the
spatial data, and \code{node1, node2} are required to fit the model without
explicitly specifying an adjacency matrix.
}
\description{
The underlying `stan` code can work without explicitly calling or building
the adjacency matrix, but to do so requires some processing, which is done
here.
}
\note{
This function is derivative from Mitzi Morris's code
\insertCite{Morris:2017,Morris:2019}{ssnet}. The "original" is
\code{mungeCARdata4stan()}, but in that function it was assumed that the
number of edges was easily calculated. Some spatial settings do not adhere
to this setup, and errors arise. This version, thus far, appears to work as
intended in less standard settings.
}
\references{
\insertRef{Morris:2017}{ssnet}

\insertRef{Morris:2019}{ssnet}
}
