#' @title Load 10x-style count matrix
#' @name load_mtx
#' @description  Loads 10x count matrices.  
#' @param mtx.folder A character to specify folder name that holds the three files (matrix.mtx - count data, features.tsv - feature data, barcodes.tsv - cell barcode data). If NULL, specific file path name has to be specified in mat.path, feat.path and bar.path. 
#' @param mat.path A character value. If mtx.folder is NULL, mat.path has to be specified for the matrix.mtx file. 
#' @param feat.path A character value. If mtx.folder is NULL, feat.path has to be specified for the features.tsv file. 
#' @param bar.path A character value. If mtx.folder is NULL, feat.path has to be specified for the barcodes.tsv file. 
#' @param is.gzip A logical. If TRUE, the input files are understood as gzipped, and handled accordingly to open.   
#' @return Returns a list consisting of gene.matrix (sparse gene count matrix), gene.features (a data.frame containing feature annotation) and barcode.features (a data.frame containing barcode annotations).
#' @export 
load_mtx <- function(mtx.folder,mat.path = NULL,feat.path = NULL,bar.path = NULL,is.gzip = FALSE)
{
  require(Matrix)
  if (!is.null(mtx.folder))
  {
    if (!is.gzip)
    {
      matrix.path = paste(mtx.folder,"/matrix.mtx",sep = "")
      features.path = paste(mtx.folder,"/features.tsv",sep = "")
      barcode.path = paste(mtx.folder,"/barcodes.tsv",sep = "")
      
      # read matrix
      mat <- readMM(file = matrix.path)
      feature.names = read.delim(features.path,header = FALSE,stringsAsFactors = FALSE)
      barcode.names = read.delim(barcode.path,header = FALSE,stringsAsFactors = FALSE)
      colnames(mat) = barcode.names$V1
      rownames(mat) = feature.names$V1
      output = list(gene.matrix = mat,gene.features = feature.names,barcode.features = barcode.names)
    }else{
      matrix.path = paste(mtx.folder,"/matrix.mtx.gz",sep = "")
      features.path = paste(mtx.folder,"/features.tsv.gz",sep = "")
      barcode.path = paste(mtx.folder,"/barcodes.tsv.gz",sep = "")
      
      # read matrix
      mat <- readMM(file = gzfile(matrix.path))
      feature.names = read.delim(gzfile(features.path),header = FALSE,stringsAsFactors = FALSE)
      barcode.names = read.delim(gzfile(barcode.path),header = FALSE,stringsAsFactors = FALSE)
      colnames(mat) = barcode.names$V1
      rownames(mat) = feature.names$V1
      output = list(gene.matrix = mat,gene.features = feature.names,barcode.features = barcode.names)
    }
  }else{
    if (!is.gzip)
    {
      # read matrix
      mat <- readMM(file = mat.path)
      feature.names = read.delim(feat.path,header = FALSE,stringsAsFactors = FALSE)
      barcode.names = read.delim(bar.path,header = FALSE,stringsAsFactors = FALSE)
      colnames(mat) = barcode.names$V1
      rownames(mat) = feature.names$V1
      output = list(gene.matrix = mat,gene.features = feature.names,barcode.features = barcode.names)
    }else{
      # read matrix
      mat <- readMM(file = gzfile(mat.path))
      feature.names = read.delim(gzfile(feat.path),header = FALSE,stringsAsFactors = FALSE)
      barcode.names = read.delim(gzfile(bar.path),header = FALSE,stringsAsFactors = FALSE)
      colnames(mat) = barcode.names$V1
      rownames(mat) = feature.names$V1
      output = list(gene.matrix = mat,gene.features = feature.names,barcode.features = barcode.names)
    }
    
  }
  
  
  return(output)
}

get_gene_qualities <- function(sce,bio.cutoff = 0,pval.cutoff = 0.2,min.cells = 10)
{
  # process data normalization
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  summary(sizeFactors(sce))
  sce = logNormCounts(sce)
  
  # evaluate 
  gene.var <- modelGeneVar(sce)
  ncell = Matrix::rowSums(counts(sce) > 1E-320,na.rm = T)
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
    dat = as.matrix(Matrix::t(logcounts(sim)))
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
