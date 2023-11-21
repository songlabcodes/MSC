##### performance metrics 
#' @title Jaccard Index Calculation 
#' @name calculate_jaccard
#' @description Calculate the Jaccard Index, as the ratio between intersection and union of two sets.  
#' @param mod.o A list containing different reference sets.  
#' @param mod.f A list containing different test sets. 
#' @return a matrix, where ith row and jth column contains the Jaccard index for ith set in mod.o and jth set in mod.f.  
#' @examples 
#' mod.o = list(set1 = LETTERS[1:10])
#' mod.f = list(set1 = LETTERS[5:20])
#' calculate_jaccard(mod.o,mod.f)
#' @export
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

#' @title Inclusion Rate Calculation 
#' @name inclusion_rate
#' @description Calculates inclusion rate.  
#' @param mod.o A list containing different reference sets.  
#' @param mod.f A list containing different test sets. 
#' @return a numeric as the overall inclusion rate.  
#' @examples 
#' mod.o = list(set1 = LETTERS[1:10])
#' mod.f = list(set1 = LETTERS[5:20])
#' inclusion_rate(mod.o,mod.f)
#' @export
inclusion_rate <- function(mod.o,mod.f)
{
  # mod.o = gold standard
  # mod.f = inferred clusters
  lst = do.call('rbind',lapply(mod.f,function(x) sapply(mod.o,function(y,z) length(intersect(y,z))/length(z),z=x)))
  ir.vec = apply(lst,1,function(x) max(x,na.rm = T))
  mlen = sapply(mod.f,length)
  sum(ir.vec*mlen)/sum(mlen)
}


#' @title Coverage Rate Calculation 
#' @name coverage_rate
#' @description Calculates coverage rate.  
#' @param mod.o A list containing different reference sets.  
#' @param mod.f A list containing different test sets. 
#' @return a numeric as the overall coverage rate.  
#' @examples 
#' mod.o = list(set1 = LETTERS[1:10])
#' mod.f = list(set1 = LETTERS[5:20])
#' coverage_rate(mod.o,mod.f)
#' @export
coverage_rate <- function(mod.o,mod.f)
{
  # mod.o = gold standard
  # mod.f = inferred clusters
  lst = do.call('rbind',lapply(mod.o,function(x) sapply(mod.f,function(y,z) length(intersect(y,z))/length(z),z=x)))
  cr.vec = apply(lst,1,function(x) max(x,na.rm = T))
  mlen = sapply(mod.o,length)
  sum(cr.vec*mlen)/sum(mlen)
}

#' @title Detection accuracy Calculation 
#' @name detection_accuracy
#' @description Calculates overall detection accuracy.  
#' @param mod.o A list containing different reference sets.  
#' @param mod.f A list containing different test sets. 
#' @return a numeric as the overall detection accuracy.  
#' @examples 
#' mod.o = list(set1 = LETTERS[1:10])
#' mod.f = list(set1 = LETTERS[5:20])
#' detection_accuracy(mod.o,mod.f)
#' @export

detection_accuracy <- function(mod.o,mod.f)
{
  jac = calculate_jaccard(mod.o,mod.f)
  ac.vec = apply(jac,1,max)
  mlen = sapply(mod.o,length)
  sum(ac.vec*mlen)/sum(mlen)
}

#### functions to do pseudo-bulk analysis
library(DESeq2)
library(Matrix.utils)
library(doParallel)


#' @title Pseudo-bulk based clusterwise differential analysis 
#' @name clusterwise_DE_analysis
#' @description Calculates differentially expressed genes (DEGs) per cell cluster by generating pseudobulk data, followed by DESeq2.  
#' @param seu Seurat object.  
#' @param modules A list object holding the list of cell names as a character vector in each element
#' @param sam.name A character value. Specify the column name in the Seurat object metadata, holding the sample assignments of individual cells.
#' @param col.name A character value. Specify the column name in the Seurat object metadata, holding case/control group assignments of individual cells.
#' @param case.id A character value. Specify the name for the case group. 
#' @param ctrl.id A character value. Specify the name for the control group. 
#' @param assay.name A character value. Specify tje assay name to use in the Seurat object.
#' @param fc.cut A numeric. Expression fold change cutoff for significance threshold. 
#' @param pval.cut A numeric. FDR cutoff for significance threshold.
#' @param n.cores A numeric. Number of cores to run the parallelized routine.
#' @return a data.frame containing the .  
#' @examples 
#' # See vignette.
#' @export

clusterwise_DE_analysis <- function(seu,modules,sam.name,
                                    case.id,ctrl.id,col.name,
                                    assay.name = "RNA",
                                    fc.cut = 2,pval.cut = 0.05,
                                    n.cores = 8)
{

  if (getDoParWorkers() == 1)
  {
    cl <- makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste0("Registered parallel cores:",getDoParWorkers(),"\n"))
  }
  
  # prepare inputs per module
  input.lst = lapply(modules,function(x,sam.name) {
    seu.sub = seu[,colnames(seu) %in%  x]
    cnt = GetAssayData(seu.sub,assay = assay.name,slot = "counts")
    meta = seu.sub@meta.data
    
    # make pseudobulk
    pseudo <- t(aggregate.Matrix(t(cnt),groupings = meta[sam.name][[1]], fun = "sum") )
    pseudo.meta = meta[match(colnames(pseudo),meta[sam.name][[1]]),];
    rownames(pseudo.meta) = colnames(pseudo)
    
    list(pb = pseudo,pmeta = pseudo.meta)
  },sam.name = sam.name  )
  pb.deg.results <- foreach (i = 1:length(input.lst),.combine = "rbind.data.frame",.packages = c("DESeq2","Matrix","Matrix.utils")) %dopar% {
    pb.results = data.frame()
    cat(paste0("processing:",names(input.lst)[i],", ",i,"/",length(input.lst),"\n"));
    # make pseudobulk
    pseudo <- as.matrix(input.lst[[i]]$pb)
    pseudo.meta = input.lst[[i]]$pmeta
    
    # remove groups with 1 count
    tbl = table(pseudo.meta[[col.name]])
    if (any(tbl <= 1))
    {
      pseudo.meta = pseudo.meta[pseudo.meta[[col.name]] %in% names(tbl)[tbl > 1],]
      pseudo = pseudo[,match(pseudo.meta[[sam.name]],colnames(pseudo))]
      pseudo.meta[[col.name]] = factor(pseudo.meta[[col.name]]) 
    }
    
    if (length(unique(pseudo.meta[[col.name]])) > 1)
    {
      mm = model.matrix(~ pseudo.meta[[col.name]] + 0)
      colnames(mm) = gsub("^(.*)\\[\\[(.*)\\]\\]","",colnames(mm))
      
      # now, get DESeq2 
      dds <- DESeqDataSetFromMatrix(data.frame(pseudo), 
                                    colData = pseudo.meta, 
                                    design = mm)
      dds <- DESeq(dds)
      
      ### get stats
      res = results(dds,contrast = list(c(case.id,ctrl.id)))
      res <- lfcShrink(dds,contrast = list(c(case.id,ctrl.id)),res=res,type = "ashr")
      res$gene.name = rownames(res)
      res$case.id = rep(case.id,nrow(res))
      res$ctrl.id = rep(ctrl.id,nrow(res))
      res$cluster.id = rep(names(modules)[i],nrow(res))
      pb.results = rbind.data.frame(pb.results,as.data.frame(res))
      rm(res)
      
    }
    return(pb.results)
  }
  
  # annotate DEGs
  pb.deg.results$DEG.ID = rep(NA,nrow(pb.deg.results))
  pb.deg.results$DEG.ID[pb.deg.results$log2FoldChange > log2(fc.cut) & pb.deg.results$padj < pval.cut] = "UP"
  pb.deg.results$DEG.ID[pb.deg.results$log2FoldChange < -log2(fc.cut) & pb.deg.results$padj < pval.cut] = "DN"
  
  return(pb.deg.results)
}