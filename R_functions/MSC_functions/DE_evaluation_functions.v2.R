library(DESeq2)
library(Matrix.utils)
library(Matrix)
library(doParallel)

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
    
    obj = subset(obj.o,cells = modules[parent.cls[i]][[1]])
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
