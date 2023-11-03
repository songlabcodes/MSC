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
  
# higher level function to take the inputs
run_msc_pipeline.v1 <- function(dat,bcnt,n.cores = 4,min.cells = 5)
{
  ### get cell network
  gf = generate_cell_network.wt.loess(mat = dat,bcnt = bcnt,
                                      dist.func = function(a, b) cor(a,b,method ="pearson"),
                                      is.decreasing = TRUE,
                                      min.size = 5,n.cores = n.cores)
  
  ## Run iterative clustering
  # identify scale parameter
  if (is.connected(gf))
  {
    alpha.res = identify_scale_parameter(g = gf,nv = 100,nr.max = 3,seed = 1234,d.func = function(x) sqrt(2*(1-x)),mode = "diameter")
    gc()
    
    iter.res = iterative_clustering.par(g.in = gf,min.size = min.cells,d.func = function(x) sqrt(2*(1-x)),alpha = alpha.res,
                                    min.width = 0.2,n.cores = n.cores)
    gc()
    iter.res$cell.network = gf;
  }else{
    cout = components(gf)
    ci = which(cout$csize > log(vcount(gf)))
    glst = lapply(ci,function(n) induced_subgraph(graph = gf,vids = names(cout$membership)[cout$membership == n]))
    
    alpha.res = mean(sapply(glst,function(gi) identify_scale_parameter(g = gi,nv = 20,nr.max = 3,seed = 1234,d.func = function(x) sqrt(2*(1-x)),mode = "diameter")))
    gc()
    
    iter.res = iterative_clustering.par(g.in = gf,min.size = min.cells,d.func = function(x) sqrt(2*(1-x)),alpha = alpha.res,
                                    min.width = 0.1,n.cores = n.cores)
  }
  iter.res$cell.network = gf;
  return(iter.res)
}


prune_hierarchy.compact <- function(iter.res)
  {
    ### prune parents for compactness gate
    parent.valid = sapply(split(iter.res$module.table$is.compact,iter.res$module.table$parent.cluster),any)
    iter.res$module.table$is.valid = iter.res$module.table$parent.cluster %in% names(parent.valid)[parent.valid] & iter.res$module.table$module.size >= 10
    sig.table = subset(iter.res$module.table,is.valid)
    
    # prune disconnected branches
    require(igraph)
    hg = graph.data.frame(sig.table[,c("cluster.name","parent.cluster")])
    cout = components(hg)
    mids = names(cout$membership)[cout$membership == cout$membership["M0"]]
    sig.table = subset(sig.table,cluster.name %in% mids)
    
    iter.res$pruned.table = sig.table
    return(iter.res)
  }
  

# function to generate single layer results
generate_clustered_data.single_layer <- function(intra.rho,inter.rho,noise.sd,cluster.size,ng)
{
  nc = sum(cluster.size)
  cluster.true = do.call('c',lapply(1:length(cluster.size),function(x) rep(x,cluster.size[x])))
  
  # create mean vector
  mu <- rep(0,nc)
  
  # create sigma matrix
  sigma = matrix(inter.rho,nrow = nc,ncol = nc)
  
  idx.f = cumsum(cluster.size)
  idx.i = c(1,idx.f[-length(idx.f)]+1)
  
  for (i in 1:length(idx.i))
  {
    sigma[idx.i[i]:idx.f[i],idx.i[i]:idx.f[i]] = intra.rho
  }
  
  diag(sigma) = 1
  
  # simulate data
  dat = mvrnorm(n=ng, mu=mu, Sigma=sigma)
  
  # add noises
  noise = matrix(rnorm(n = nc * ng,mean = 0,sd = noise.sd),nrow = ng,ncol = nc)
  dat = dat + noise
  
  # add dimnames
  colnames(dat) = paste0("cell",1:ncol(dat))
  rownames(dat) = paste0("gene",1:nrow(dat))
  names(cluster.true) = colnames(dat)
  
  # generate binary matrix to illuminate presence of valid data 
  bcnt = matrix(1,nrow = nrow(dat),ncol = ncol(dat))
  dimnames(bcnt) = dimnames(dat)
  
  return(list(data = dat,binary.data = bcnt,cluster.true = cluster.true))
  
}

generate_clustered_data.multi_layer <- function(cluster.size,merge.instr,rho.vec,ng,noise.sd)
{
  nc = sum(cluster.size)
  
  # create mean vector
  mu <- rep(0,nc)
  
  # create sigma matrix
  sigma = matrix(NA,nrow = nc,ncol = nc)
  cluster.o = do.call('c',lapply(1:length(cluster.size),function(x) rep(x,cluster.size[x])))
  cluster.matrix= matrix(0,nrow = nc,ncol = length(merge.instr))
  for (i in 1:length(merge.instr))
  {
    cluster.merge = rep(NA,nc)
    cls.o = unique(cluster.o)
    for (j in 1:length(merge.instr[[i]])) cluster.merge[cluster.o == cls.o[j]] = merge.instr[[i]][j]
    cluster.matrix[,i] = cluster.merge;
    
    # put in rho
    cls.merge = unique(cluster.merge)
    for (j in 1:length(cls.merge))
    {
      ii = which(cluster.merge == cls.merge[j])
      m = sigma[ii,ii]
      m[is.na(m)] = rho.vec[i]
      sigma[ii,ii] = m
      rm(m)
    }
    
    cluster.o = cluster.merge
    rm(cluster.merge)
  }
  
  # set diagonal 
  diag(sigma) = 1
  
  # simulate data
  dat = mvrnorm(n=ng, mu=mu, Sigma=sigma)
  
  # add noises
  noise = matrix(rnorm(n = nc * ng,mean = 0,sd = noise.sd),nrow = ng,ncol = nc)
  dat = dat + noise
  
  # add dimnames
  colnames(dat) = paste0("cell",1:ncol(dat))
  rownames(dat) = paste0("gene",1:nrow(dat))
  rownames(cluster.matrix) = colnames(dat)
  
  # generate binary matrix to illuminate presence of valid data 
  bcnt = matrix(1,nrow = nrow(dat),ncol = ncol(dat))
  dimnames(bcnt) = dimnames(dat)
  
  return(list(data = dat,binary.data = bcnt,cluster.matrix = cluster.matrix))
}

calculate_jaccard <- function(mod.o,mod.f)
{
  jac = matrix(0,nrow = length(mod.o),ncol = length(mod.f));
  rownames(jac) = names(mod.o)
  colnames(jac) = names(mod.f)
  for (i in 1:length(mod.o))
  {
    len = length(mod.o[[i]])
    xo = mod.o[[i]]
    jac[i,] = sapply(mod.f,function(x) {length(intersect(x,xo))/length(union(x,xo))})
  }
  return(jac)
}
