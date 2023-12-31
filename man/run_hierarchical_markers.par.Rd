% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DE_evaluation_functions.v3.R
\name{run_hierarchical_markers.par}
\alias{run_hierarchical_markers.par}
\title{Perform hierarchical cluster marker detection}
\usage{
run_hierarchical_markers.par(
  obj.o,
  mtbl,
  modules,
  grp.var,
  cls.var = "msc.cluster",
  find.conserved = FALSE,
  assay.name = "SCT",
  recorrect.SCT = FALSE,
  DE.method = "MAST",
  lfc.cutoff = log2(1.2),
  min.size = 50,
  n.cores = 6
)
}
\arguments{
\item{obj.o}{A Seurat object.}

\item{mtbl}{A data.frame for MSC module table, capturing child and parent cluster (mtbl$parent.cluster) relationship.}

\item{modules}{A list object holding the list of cell names as a cell cluster in each.}

\item{grp.var}{Column name in Seurat meta data, designating sample assignments.}

\item{cls.var}{Column name in Seurat meta data, designating child cluster assignments. "msc.cluster" is the default, and this column is created within this function to evaluate each split per parent cluster.}

\item{find.conserved}{A logical. If TRUE, conserved markers across samples specific in "grp.var" column are calculated.}

\item{assay.name}{Seurat assay name to perform the marker analysis.}

\item{recorrect.SCT}{A logical. If TRUE, perform PrepSCTFindMarkers() to re-scale the data across different samples specific in "grp.var".}

\item{DE.method}{Differential expression test method, as specified in \code{\link[=FindMarkers]{FindMarkers()}} in Seurat package.}

\item{lfc.cutoff}{log2 fold change cutoff. Default is log2(1.2).}

\item{min.size}{An integer to specify the minumum cluster size to compute the markers.}

\item{n.cores}{An integer for the number of cores for parallelization}
}
\value{
Returns a named list. Each element contains a data.frame for list of markers per the named cluster.
}
\description{
Wrapper of FindMarkers() function in Seurat to perform marker detection comparing child clusters, derived from their parent cluster.
}
\examples{
\donttest{
data(simMix1)

library(scran)
library(igraph)
library(doParallel)
qc.out = process_data(simMix1,do.impute = TRUE)

# add reduced dim
reducedDim(qc.out$sce,"PCA.imputed") = calculatePCA(x = assay(qc.out$sce,"ALRA.logcounts"),subset_row = rownames(subset(qc.out$gene.data.alra,p.value < 0.05)))
reducedDim(qc.out$sce,"UMAP.imputed") = calculateUMAP(x = qc.out$sce,dimred = "PCA.imputed")
reducedDim(qc.out$sce,"tSNE.imputed") = calculateTSNE(x = qc.out$sce,dimred = "PCA.imputed")

plotReducedDim(qc.out$sce,dimred = "tSNE.imputed",colour_by = "phenoid")

dat = subset(qc.out$gene.data.alra,bio > 0 & p.value < 0.05)
bcnt.alra = counts(qc.out$sce)[rownames(dat),]
m.alra = as.matrix(assay(qc.out$sce,"ALRA.logcounts")[rownames(dat),])
bcnt.alra[bcnt.alra > 1E-320] = 1
msc.res = msc_workflow(dat = assay(qc.out$sce,"ALRA.logcounts"),bcnt = bcnt.alra,n.cores = 4)
msc.res = prune_hierarchy.compact(msc.res)

## conversion to Seurat object
library(Seurat)
seu = CreateSeuratObject(counts = counts(qc.out$sce),project = "simMix1",meta.data = as.data.frame(colData(qc.out$sce)))
seu@assays$RNA@data = assay(qc.out$sce,"ALRA.logcounts")

msc.marker.list = run_hierarchical_markers.par(obj.o = seu,mtbl = msc.res$pruned.table,
modules = msc.res$modules,grp.var = "orig.ident",cls.var = "msc.cluster",find.conserved = FALSE,assay.name = "RNA",recorrect.SCT = FALSE,
DE.method = "wilcox",n.cores = 4)
}
}
