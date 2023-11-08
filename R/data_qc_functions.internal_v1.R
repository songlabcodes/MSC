get_gene_qualities <- function(sce,bio.cutoff = 0,pval.cutoff = 0.2,min.cells = 10)
{
  # process data normalization
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  summary(sizeFactors(sce))
  sce = logNormCounts(sce)
  
  # evaluate 
  gene.var <- modelGeneVar(sce)
  ncell = rowSums(counts(sce) > 1E-320,na.rm = T)
  gene.var$cell.express = ncell;
  gene.var$bio.filtered = rownames(gene.var) %in% rownames(subset(gene.var,bio > bio.cutoff & p.value < pval.cutoff & cell.express >= min.cells ))
  
  #### Now, get variable genes to use 
  sig.gene.var = subset(gene.var,bio.filtered)
  dgenes = rownames(sig.gene.var)
  
  cat(paste0("- number of genes filtered:",length(dgenes),"\n"))
  return(gene.var)
}

#' @title Basic data processing with scra
#' @name process_data
#' @description  Use scran workflow to get variable genes, then use ALRA for imputation. 
#' @param sim A SingleCellExperiment object with counts in it.
#' @param do.impute A logical. If TRUE, ALRA imputation is performed. Note that if data is big, the imputation takes time.   
#' @return Returns a list with processed SingleCellExperiment object, data.frame for gene statistics, then a vector for listing variable genes.
#' @examples 
#' \donttest{
#' data(simMix1)
#' library(scran)
#' library(scater)
#' qc.out = process_data(simMix1,do.impute = TRUE)
#' 
#' # add reduced dim
#' reducedDim(qc.out$sce,"PCA.imputed") = calculatePCA(x = assay(qc.out$sce,"ALRA.logcounts"),subset_row = rownames(subset(qc.out$gene.data.alra,p.value < 0.05)))
#' reducedDim(qc.out$sce,"UMAP.imputed") = calculateUMAP(x = qc.out$sce,dimred = "PCA.imputed")
#' reducedDim(qc.out$sce,"tSNE.imputed") = calculateTSNE(x = qc.out$sce,dimred = "PCA.imputed")
#' 
#' plotReducedDim(qc.out$sce,dimred = "tSNE.imputed",colour_by = "phenoid")
#' }
#' @export
process_data <- function(sim,do.impute = FALSE)
{
  # process data normalization
  clusters <- quickCluster(sim)
  sce <- computeSumFactors(sim, clusters=clusters)
  summary(sizeFactors(sim))
  sim = logNormCounts(sim)
  
  gene.var = get_gene_qualities(sce = sim,bio.cutoff = 0,pval.cutoff = 0.2,min.cells = 10)
  dgenes = getTopHVGs(stats = gene.var,n = 3000)
  
  # denoise PCs
  sim <- fixedPCA(sim, subset.row=dgenes)
  sim <- denoisePCA(sim, gene.var, subset.row=dgenes)
  ncol(reducedDim(sim, "PCA"))
  
  output <- getClusteredPCs(reducedDim(sim))
  npcs <- metadata(output)$chosen
  reducedDim(sim, "PCAsub") <- reducedDim(sim, "PCA")[,1:npcs,drop=FALSE]
  
  # update tSNE
  library(scater)
  sim <- runTSNE(sim, dimred="PCAsub")
  sim <- runUMAP(sim, dimred="PCAsub")
  
  # run imputation: ALRA
  if (do.impute)
  {
    dat = as.matrix(t(logcounts(sim)))
    alra.res = alra(A_norm= dat);rm(dat)
    m = t(alra.res[[3]])
    dimnames(m) = dimnames(logcounts(sim))
    alra.var = modelGeneVar(x = m)
    assay(sim,"ALRA.logcounts") = m
    output = list(sce = sim,gene.data = gene.var,variable.genes = dgenes,gene.data.alra = alra.var)
    
  }else{
    output = list(sce = sim,gene.data = gene.var,variable.genes = dgenes)
    
  }
  
  
  
  return(output)
  #plot_grid(plotUMAP(sim,colour_by = "Group"),plotTSNE(sim,colour_by = "Group"))
}
