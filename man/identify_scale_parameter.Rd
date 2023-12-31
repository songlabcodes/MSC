% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_functions.v3.4.R
\docType{package}
\name{identify_scale_parameter}
\alias{identify_scale_parameter}
\title{Calculate the alpha value for compactness}
\arguments{
\item{g}{An igraph object containing the network}

\item{nv}{An integer. Number of random nodes to sample for neighborhood search.}

\item{nr.max}{Neighborhood layer to expand.}

\item{d.func}{a function to convert the edge weights into edge function.}

\item{mode}{A character value from either of "diameter" or "mean". Decide to estimate the compactness based on network diameter (i.e. the longest of the pairwise shortest path distance), or mean (i.e. mean of all pairwise shortest path distance).}
}
\value{
Returns a numeric value for the alpha value.
}
\description{
Estimates the compactness parameter alpha via subnetwork sampling.
}
\examples{
data(pbmc_8k_msc_results)
identify_scale_parameter(pbmc_8k_msc_results$cell.network,mode = "diameter")
}
