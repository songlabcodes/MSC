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
  
  process_data <- function(sim)
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
    dat = as.matrix(t(logcounts(sim)))
    alra.res = alra(A_norm= dat);rm(dat)
    m = t(alra.res[[3]])
    dimnames(m) = dimnames(logcounts(sim))
    alra.var = modelGeneVar(x = m)
    assay(sim,"ALRA.logcounts") = m
    
    
    output = list(sce = sim,gene.data = gene.var,variable.genes = dgenes,gene.data.alra = alra.var)
    return(output)
    #plot_grid(plotUMAP(sim,colour_by = "Group"),plotTSNE(sim,colour_by = "Group"))
  }
  
  get_snn_clusters <- function(seu,dimred,dims,resol.vec = seq(0.4,1.2,0.1))
  {
    seu <- FindNeighbors(seu,reduction = dimred, dims = dims)
    seu <- FindClusters(seu, resolution = resol.vec)
    seu.cls.col = colnames(seu@meta.data)[grep("snn_res",colnames(seu@meta.data))]
    
    snn.data = seu@meta.data[,seu.cls.col]
    rownames(snn.data) = colnames(seu)
    colnames(snn.data) = gsub("^(.*)_snn_res","seurat_snn_res",colnames(snn.data))
    snn.mods = lapply(snn.data,function(x) split(rownames(snn.data),factor(x)))
    
    for (i in 1:length(seu.cls.col))
    {
      names(snn.mods[[i]]) = paste0(names(snn.mods)[i],"_C",names(snn.mods[[i]]))
    }
    names(snn.mods) = NULL
    snn.mods = do.call('c',snn.mods)
    
    output = list(modules = snn.mods,cluster.matrix = snn.data)
    
    return(output)
  }
  
  
  
  
