% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lib_iKernel.R
\name{get_voronoi_feature}
\alias{get_voronoi_feature}
\alias{add_new_point_iKernel}
\title{The function to get feature representation in RKHS based on Voronoi diagram for WHOLE dataset}
\usage{
get_voronoi_feature(
  psi = 40,
  t = 350,
  data,
  talkative = FALSE,
  Matrix_Voronoi = NULL
)

add_new_point_iKernel(data, d1, Matrix_Voronoi, dissim, t, psi, nr)
}
\arguments{
\item{psi}{Integer number related to the size of each Voronoi diagram}

\item{t}{Integer number of trees in Isolation Kernel or dimension of RKHS}

\item{data}{dataset of points, rows - points, columns - dimensions of a point}

\item{talkative}{logical. If TRUE then print messages, FALSE for the silent execution}

\item{Matrix_Voronoi}{Matrix of Voronoi diagrams, if it is NULL then the function will calculate Matrix_Voronoi}

\item{d1}{Data point - usually it is an observation data point}

\item{dissim}{Matrix of dissimilarity or distances between all points.}

\item{nr}{Integer number of rows in matrix of distances (dissim) and also the size of dataset}
}
\value{
Feature representation in RKHS based on Voronoi diagram for WHOLE dataset

RKHS mapping for a new point based on Isolation Kernel mapping
}
\description{
The function to get feature representation in RKHS based on Voronoi diagram for WHOLE dataset
}
\section{Functions}{
\itemize{
\item \code{add_new_point_iKernel()}: The function to get RKHS mapping based on Isolation Kernel for a new point

}}
\examples{
NULL
NULL
}
\keyword{internal}
