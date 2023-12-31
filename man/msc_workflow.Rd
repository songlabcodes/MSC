% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_functions.v3.4.R
\docType{package}
\name{msc_workflow}
\alias{msc_workflow}
\title{multi-scale clustering workflow}
\usage{
msc_workflow(
  dat,
  bcnt,
  n.cores = 4,
  min.cells = 5,
  sim.func = function(a, b) cor(a, b, method = "pearson"),
  is.decreasing = TRUE,
  d.func = function(x) sqrt(2 * (1 - x)),
  valid.gate = "both"
)
}
\arguments{
\item{dat}{A log-normalized single-cell expression matrix (rows: genes, columns: cell).}

\item{bcnt}{A binary matrix of genes (row) by cells (column). If expressed, \code{bcnt[i,j] = 1}, otherwise 0.}

\item{n.cores}{An integer for the number of cores for parallelization}

\item{min.cells}{An integer specifying the minumum cell cluster size}

\item{sim.func}{A function to calculate cell-cell similarity. Default is the pairwise Pearson's correlation}

\item{is.decreasing}{A logical. TRUE if dist.func is similarity function to show greater value corresponds to greater similarity. FALSE if dist.func is dissimilarity.}

\item{d.func}{A function to convert the pairwise similarity to pairwise distance.}

\item{valid.gate}{Filtering method to identify significant clusters. "density" compares intra-cluster link density between parent and child clusters. "compact" compared compactness with the parent clusters. "both" uses "density" and "compact" filters simultaneously. Default is "both".}
}
\value{
Returns a list containing clustering results with optimal.modules (the final results from adapt split), resolution.tested (clustering resolution sweeped), resolution.used (the final value of clustering resolution), and module.table containing the compactness and intra-cluster connectivity stats.
}
\description{
Wrapper function to perform LEN and top-down iterative clustering altogether.
}
\examples{
\donttest{
data(simMix1)
library(scater)
library(scran)
library(igraph)
library(doParallel)
qc.out = process_data(simMix1,do.impute = TRUE)

# add reduced dim
reducedDim(qc.out$sce,"PCA.imputed") = calculatePCA(x = assay(qc.out$sce,"ALRA.logcounts"),subset_row = rownames(subset(qc.out$gene.data.alra,p.value < 0.05)))
reducedDim(qc.out$sce,"UMAP.imputed") = calculateUMAP(x = qc.out$sce,dimred = "PCA.imputed")
reducedDim(qc.out$sce,"tSNE.imputed") = calculateTSNE(x = qc.out$sce,dimred = "PCA.imputed")

plotReducedDim(qc.out$sce,dimred = "tSNE.imputed",colour_by = "phenoid")
# for imputed data for correlation computation 
dat = subset(qc.out$gene.data.alra,bio > 0 & p.value < 0.05)
bcnt.alra = counts(qc.out$sce)[rownames(dat),]
m.alra = as.matrix(assay(qc.out$sce,"ALRA.logcounts")[rownames(dat),])
bcnt.alra[bcnt.alra > 1E-320] = 1
msc.res = msc_workflow(dat = m.alra,bcnt = bcnt.alra,n.cores = 4)
}
}
