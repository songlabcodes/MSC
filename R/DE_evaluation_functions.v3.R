
### marker functions
annotate_conserved_markers <- function(marker.list,pval.cutoff = 0.05,fc.cutoff = 1.2,pct.1.cutoff = 0.25,pct.diff.cutoff = 0.1)
{
  require(ACAT)
  marker.list = marker.list[!sapply(marker.list,is.null)]
  
  if (length(marker.list) == 0) stop("All marker results are NULL.")
  marker.sig <- lapply(marker.list,function(x,pval,lfc,pct.1,pct.diff,nr) {
    # get the adjusted p-value matrix
    adj.pval <- as.matrix(x[,grep("p_val_adj$",colnames(x))])
    adj.pval.sig <- apply(adj.pval,1,function(p,q) all(p < q),q = pval)
    
    # get summarized p-value
    p.sum = ACAT(t(as.matrix(x[,grep("p_val$",colnames(x))])))
    #adjp.sum = p.sum * nr;adjp.sum[adjp.sum > 1] = 1
    adjp.sum = p.adjust(p.sum,"BH");
    
    # get logfc matrix
    lfcm = as.matrix(x[,grep("avg_log2FC$",colnames(x))])
    lfcm.sig = apply(lfcm,1,function(p,q) all(p > q),q = lfc)
    lfc.avg = rowMeans(lfcm,na.rm = TRUE)
    
    # get % expressed
    pct.1.m = as.matrix(x[,grep("pct.1$",colnames(x))])
    pct.1.sig = apply(pct.1.m,1,function(p,q) all(p > q),q = pct.1)
    
    # get % differently expressed 
    cohort = gsub("_pct.1","",colnames(x)[grep("pct.1$",colnames(x))])
    pct.diff.m = do.call('cbind',lapply(cohort,function(p,tbl) {m = as.matrix(tbl[,match(paste0(p,c("_pct.1","_pct.2")),colnames(tbl))]);m[,1] - m[,2]},tbl = x))
    colnames(pct.diff.m) = cohort
    pct.diff.sig = apply(pct.diff.m,1,function(p,q) all(p > q),q = pct.diff)
    
    out = x;
    out$p_val.sum = p.sum
    out$adj_p_val.sum = adjp.sum
    out$adj_p_val.sig = adj.pval.sig
    out$lfc.sig = lfcm.sig
    out$lfc.avg = lfc.avg
    out$pct.1.sig = pct.1.sig
    out$pct.diff.sig = pct.diff.sig
    out$gene.symbol = rownames(out)
    out$overall.sig = ( adjp.sum < pval.cutoff & lfcm.sig & pct.1.sig & pct.diff.sig)
    
    return(out)
  },pval = pval.cutoff,lfc = log2(fc.cutoff),pct.1 = pct.1.cutoff,pct.diff = pct.diff.cutoff)
  return(marker.sig)
}

run_conserved_markers.par <- function(obj,find.conserved = TRUE,
                                      grp.var,cls.var,DE.method = "MAST",
                                      lfc.cutoff = log2(1.2),
                                      assay.name = "SCT",recorrect.SCT = FALSE,n.cores = 4)
{
  require(doParallel)
  n.reg = getDoParWorkers()
  
  if (n.reg == 1 & n.cores > 1)
  {
    cl = makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Registered cores:",getDoParWorkers(),"\n"))
  }
  
  Idents(obj) = factor(obj@meta.data[[which(colnames(obj@meta.data) == cls.var)]])
  cls = levels(Idents(obj))
  
  if (assay.name == "SCT" & recorrect.SCT) obj = PrepSCTFindMarkers(obj)
  
  if (find.conserved)
  {
    recorrect_umi = FALSE
    if (recorrect.SCT & assay.name == "SCT") recorrect_umi = FALSE
    
    job.idx = do.call('c',lapply(1:n.cores,function(x) rep(x,ceiling(length(cls)/n.cores))))[1:length(cls)]
    cls.idx = split(cls,factor(job.idx))
    
    marker.list <- foreach(k.idx = cls.idx,.export = c("grp.var","recorrect_umi","lfc.cutoff"),.packages = c("Seurat","ACAT"),.combine = 'c') %dopar% {
      mlst = vector("list",length(k.idx))
      names(mlst) = k.idx
      for (k in 1:length(k.idx))
      {
        out = try(FindConservedMarkers(object = obj,assay = assay.name,ident.1 = k.idx[k],grouping.var = grp.var,
                                       test.use = DE.method,recorrect_umi = recorrect_umi,
                                       logfc.threshold = lfc.cutoff,min.cells.group = 10));
        if(!inherits(out, 'try-error'))  
        {
          out$gene.name = rownames(out)
          out$cluster.id = rep(k.idx[k],nrow(out))
          
        }else{
          out = NULL
        }
        
        mlst[[k]] = out
        
      }
      
      return(out)
    }
    marker.list.annot = annotate_conserved_markers(marker.list,pval.cutoff = 0.05,fc.cutoff = 2^lfc.cutoff,pct.1.cutoff = 0.25,pct.diff.cutoff = 0.1)
    
  }else{
    recorrect_umi = FALSE
    if (recorrect.SCT & assay.name == "SCT") recorrect_umi = FALSE
    
    job.idx = do.call('c',lapply(1:n.cores,function(x) rep(x,ceiling(length(cls)/n.cores))))[1:length(cls)]
    cls.idx = split(cls,factor(job.idx))
    
    marker.list.annot <- foreach(k.idx = cls.idx,.export = c("grp.var","recorrect_umi","lfc.cutoff"),.packages = c("Seurat","ACAT"),.combine = 'c') %dopar% {
      mlst = vector("list",length(k.idx))
      names(mlst) = k.idx
      for (k in 1:length(k.idx))
      {
        out = try(FindMarkers(object = obj,assay = assay.name,ident.1 = k.idx[k],test.use = DE.method,
                              recorrect_umi = recorrect_umi,
                              logfc.threshold = lfc.cutoff,min.cells.group = 10));
        
        if(!inherits(out, 'try-error'))  
        {
          out$gene.name = rownames(out)
          out$cluster.id = rep(k.idx[k],nrow(out))
        }else{
          out = NULL
        }
        mlst[[k]] = out
        rm(out)
      }
      
      
      return(mlst)
    }
    
  }
  return(marker.list.annot)
}

#' @title Perform hierarchical cluster marker detection 
#' @name run_hierarchical_markers.par
#' @description  Wrapper of FindMarkers() function in Seurat to perform marker detection comparing child clusters, derived from their parent cluster. 
#' @param obj.o A Seurat object.
#' @param mtbl A data.frame for MSC module table, capturing child and parent cluster (mtbl$parent.cluster) relationship.
#' @param modules A list object holding the list of cell names as a cell cluster in each.   
#' @param grp.var Column name in Seurat meta data, designating sample assignments.
#' @param cls.var Column name in Seurat meta data, designating child cluster assignments. "msc.cluster" is the default, and this column is created within this function to evaluate each split per parent cluster. 
#' @param find.conserved A logical. If TRUE, conserved markers across samples specific in "grp.var" column are calculated. 
#' @param DE.method Differential expression test method, as specified in [FindMarkers()] in Seurat package.
#' @param lfc.cutoff log2 fold change cutoff. Default is log2(1.2).
#' @param assay.name Seurat assay name to perform the marker analysis.
#' @param min.size An integer to specify the minumum cluster size to compute the markers. 
#' @param recorrect.SCT A logical. If TRUE, perform PrepSCTFindMarkers() to re-scale the data across different samples specific in "grp.var". 
#' @param n.cores An integer for the number of cores for parallelization 
#' @return Returns a named list. Each element contains a data.frame for list of markers per the named cluster.
#' @examples
#' \donttest{
#' data(simMix1)
#' 
#' library(scran)
#' library(igraph)
#' library(doParallel)
#' qc.out = process_data(simMix1,do.impute = TRUE)
#' 
#' # add reduced dim
#' reducedDim(qc.out$sce,"PCA.imputed") = calculatePCA(x = assay(qc.out$sce,"ALRA.logcounts"),subset_row = rownames(subset(qc.out$gene.data.alra,p.value < 0.05)))
#' reducedDim(qc.out$sce,"UMAP.imputed") = calculateUMAP(x = qc.out$sce,dimred = "PCA.imputed")
#' reducedDim(qc.out$sce,"tSNE.imputed") = calculateTSNE(x = qc.out$sce,dimred = "PCA.imputed")
#' 
#' plotReducedDim(qc.out$sce,dimred = "tSNE.imputed",colour_by = "phenoid")
#' 
#' dat = subset(qc.out$gene.data.alra,bio > 0 & p.value < 0.05)
#' bcnt.alra = counts(qc.out$sce)[rownames(dat),]
#' m.alra = as.matrix(assay(qc.out$sce,"ALRA.logcounts")[rownames(dat),])
#' bcnt.alra[bcnt.alra > 1E-320] = 1
#' msc.res = msc_workflow(dat = assay(qc.out$sce,"ALRA.logcounts"),bcnt = bcnt.alra,n.cores = 4)
#' msc.res = prune_hierarchy.compact(msc.res)
#' 
#' ## conversion to Seurat object
#' library(Seurat)
#' seu = CreateSeuratObject(counts = counts(qc.out$sce),project = "simMix1",meta.data = as.data.frame(colData(qc.out$sce)))
#' seu@assays$RNA@data = assay(qc.out$sce,"ALRA.logcounts")
#' 
#' msc.marker.list = run_hierarchical_markers.par(obj.o = seu,mtbl = msc.res$pruned.table,
#' modules = msc.res$modules,grp.var = "orig.ident",cls.var = "msc.cluster",find.conserved = FALSE,assay.name = "RNA",recorrect.SCT = FALSE,
#' DE.method = "wilcox",n.cores = 4)
#' }
#' @export
run_hierarchical_markers.par <- function(obj.o,mtbl,modules,grp.var,cls.var = "msc.cluster",
                                         find.conserved = FALSE,assay.name = "SCT",recorrect.SCT = FALSE,
                                         DE.method = "MAST",lfc.cutoff = log2(1.2),
                                         min.size = 50,
                                         n.cores = 6)
{
  
  require(doParallel)
  n.reg = getDoParWorkers()
  
  if (n.reg == 1 & n.cores > 1)
  {
    cl = makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Registered cores:",getDoParWorkers(),"\n"))
  }
  
  ### generate results per root.
  marker.results = list()
  
  # get parent with valid 
  parent.cls = unique(mtbl$parent.cluster)
  parent.cls = parent.cls[sapply(parent.cls,function(x) nrow(subset(mtbl,parent.cluster == x & module.size >= min.size)) > 1)]
  
  for (i in 1:length(parent.cls))
  {
    cat(paste0("processing parent cluster:",parent.cls[i],"\n"));
    
    #
    selected <- WhichCells(obj.o,cells = modules[parent.cls[i]][[1]])
    obj = subset(obj.o,cells = selected)
    print(dim(obj))
    
    # get subclusters
    mid.f = modules[subset(mtbl,parent.cluster == parent.cls[i])$cluster.name]
    cat(paste0("- child clusters are:",paste(names(mid.f),collapse = ","),"\n"))
    
    vec = rep(NA,ncol(obj));names(vec) = colnames(obj)
    for (m in 1:length(mid.f)) vec[names(vec) %in% mid.f[[m]]] = names(mid.f)[m]
    obj@meta.data[cls.var] = factor(vec)
    
    # run cluster routine
    sig.marker.list = run_conserved_markers.par(obj = obj,find.conserved = find.conserved,
                                                grp.var = grp.var,cls.var = cls.var,
                                                lfc.cutoff = lfc.cutoff,
                                                assay.name = assay.name,recorrect.SCT = recorrect.SCT,
                                                DE.method = DE.method,
                                                n.cores = n.cores)
    marker.results = c(marker.results,sig.marker.list)
    rm(obj,mid.f,vec,sig.marker.list)
    gc()
  }
  return(marker.results)
}

get_dotplot_pdata <- function(genes,cnt,mat,grp,cls)
{
  # genes =list of gene names for plotting
  # cnt = count matrix
  # mat = normalized expression matrix
  # grp = character vector specifying group assignments for cells (e.g. case / control) 
  # cls = character vector specifying cluster assignments for cells
  gv = unique(grp)
  cv = unique(cls)
  
  nc <- grp.vec <- cls.vec <- pct.expr <- mean.expr <- gene.name <- c() 
  for (i in 1:length(genes))
  {
    cnt.g = as.vector(cnt[match(genes[i],rownames(cnt)),]);names(cnt.g) = colnames(cnt)
    mat.g = as.vector(mat[match(genes[i],rownames(mat)),]);names(mat.g) = colnames(mat)
    for (n in 1:length(gv))
    {
      for (m in 1:length(cv))
      {
        ii = which(grp == gv[n] & cls == cv[m])
        if (length(ii) > 0)
        {
          nc = c(nc,length(ii))
          gene.name <- c(gene.name,genes[i])
          grp.vec <- c(grp.vec,gv[n])
          cls.vec <- c(cls.vec,cv[m])
          pct.expr <- c(pct.expr,sum(cnt.g[ii] > 1E-320,na.rm = T)/length(ii) * 100)
          mean.expr <- c(mean.expr,mean(mat.g[ii]))
        }
      }
    }
  }
  out = data.frame(gene.name = gene.name,group.name = grp.vec,cluster.name = cls.vec,pct.express = pct.expr,mean.express = mean.expr,num.cell = nc)
  return(out)
}
