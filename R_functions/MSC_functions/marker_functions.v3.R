# get coordinates for centered labels 
get_module_label_coord <- function(coord,mods)
{
	do.call('rbind',lapply(mods,function(x,y) apply(y[match(x,rownames(y)),],2,function(q) median(q,na.rm = T)),y = coord))
}
create_pseudobulk_per_cluster <- function(seur,mods,group.var,
                                          assay.name = "RNA")
{
  require(Matrix.utils)
  vec.c = rep(NA,length(mods))
  for (m in 1:length(mods)) vec.c[colnames(seur) %in% mods[[m]]] = names(mods)[m]
  seur@meta.data$msc_clusters = vec.c
  
  groups = seur@meta.data[,c("msc_clusters",group.var)]
  pb <- aggregate.Matrix(t(GetAssayData(seur[[assay.name]],slot = "counts")), 
                         groupings = groups, fun = "sum") 
  
  pb = t(pb)
  pb = pb[,which(colnames(pb) != "NA")]
  #pb = pb + 1
  
  # remove unexpressed genes
  ii = which(rowSums(pb < 1E-320,na.rm = TRUE) < ceiling(ncol(pb)*0.8))
  pb = pb[ii,]
  print(dim(pb))
  
  # create sample meta data
  sample.meta = data.frame(id = colnames(pb),cluster.id = factor(gsub("_(.*)$","",colnames(pb))),sample.id = factor(gsub("^(.*)_","",colnames(pb))))
  
  return(list(pseudobulk = pb,meta = sample.meta))
}


find_pseudobulk_markers <- function(count.mat,meta,
                                    seur,modules,
                                    assay.name = "RNA",
                                    pval.cut = 0.05,fc.cut = 1.2,
                                    do.par = FALSE,n.cores = 8)
{
  # count.mat = count matrix in sparse format, meta = data frame with columns "cluster.id" and "sample.id". 
  dds <- tryCatch(DESeqDataSetFromMatrix(count.mat, colData = meta, design = ~ 0  + cluster.id + sample.id),error = function(e) return(NA))
  
  # Run DESeq2 differential expression analysis
  if (do.par) 
  {
    gc()
    mc.param = MulticoreParam(workers = n.cores)
    register(mc.param)
    dds.o <- tryCatch(DESeq(dds,parallel = TRUE,BPPARAM = mc.param),error = function(e) return(NA))
    if (any(is.na(dds.o)) )
    {
      dds = DESeq(dds,parallel = TRUE,fitType = "mean",BPPARAM = mc.param)
    }else{
      dds = dds.o
    }
    rm(dds.o)
  }
  
  if (!do.par) 
  {
    
    dds.o <- tryCatch(DESeq(dds),error = function(e) return(NA))
    if (any(is.na(dds.o)) )
    {
      dds = DESeq(dds,fitType = "mean")
    }else{
      dds = dds.o
    }
    rm(dds.o)
  }
  gc()
  # extract results for genotype comparison
  cls.var = resultsNames(dds);
  cls.var = cls.var[grep("^cluster.id",cls.var)]
  print(gsub("cluster.id","",cls.var))
  
  if (!do.par)
  {
    cluster.res = vector("list",length(cls.var))
    names(cluster.res) = gsub("cluster.id","",cls.var)
    for (j in 1:length(cls.var))
    {
      case = cls.var[j];ctrl = cls.var[-j]
      res = results(dds,contrast = list(case,ctrl),listValues = c(1,-1/length(ctrl)))
      res <- lfcShrink(dds,contrast = list(case,ctrl),res=res,type = "ashr")
      
      # call results
      res$gene.name = rownames(res)
      res$case.group = rep(gsub("cluster.id","",case),nrow(res))
      res$ctrl.group = rep(paste(gsub("cluster.id","",ctrl),collapse= ","),nrow(res))
      
      # call DEG
      deg = rep(NA,nrow(res))
      deg[res$log2FoldChange > log2(fc.cut) & res$padj < pval.cut] = "UP"
      deg[res$log2FoldChange < -log2(fc.cut) & res$padj < pval.cut] = "DN"
      res$DEG.ID = deg;
      
      # add percent expressed
      m = GetAssayData(seur,assay= assay.name,slot = "counts")
      cc = which(colnames(m) %in% modules[gsub("^cluster.id","",cls.var[j])][[1]])
      cc.not = which(!(colnames(m) %in% modules[gsub("^cluster.id","",cls.var[j])][[1]]))
      m.in = m[match(res$gene.name,rownames(m)),cc]
      m.not = m[match(res$gene.name,rownames(m)),cc.not]
      
      res$pct.case = rowSums(m.in > 1E-320,na.rm = TRUE)/ncol(m)
      res$pct.ctrl = rowSums(m.not > 1E-320,na.rm = TRUE)/ncol(m)
      #res$parent.module = rep(vr[i],nrow(res))
      cluster.res[[j]] = res
      
      rm(m,cc,cc.not,m.in,m.not,deg,case,ctrl)
    }
  }else{
    require(foreach)
    require(iterators)
    require(doParallel)
    if (getDoParWorkers() < n.cores) 
    {
      # remove pre-registered nodes to make sure it runs clean.
      if (getDoParWorkers() > 1)
      {
        env <- foreach:::.foreachGlobals
        rm(list = ls(name = env),pos = env)
      }
      cl <- parallel::makeCluster(n.cores)
      doParallel::registerDoParallel(cl)
    }
    cluster.res = foreach(j=1:length(cls.var),.packages = c("DESeq2","Seurat") ) %dopar% {
      case = cls.var[j];ctrl = cls.var[-j]
      res = results(dds,contrast = list(case,ctrl),listValues = c(1,-1/length(ctrl)))
      res <- lfcShrink(dds,contrast = list(case,ctrl),res=res,type = "ashr")
      
      # call results
      res$gene.name = rownames(res)
      res$case.group = rep(gsub("cluster.id","",case),nrow(res))
      res$ctrl.group = rep(paste(gsub("cluster.id","",ctrl),collapse= ","),nrow(res))
      
      # call DEG
      deg = rep(NA,nrow(res))
      deg[res$log2FoldChange > log2(fc.cut) & res$padj < pval.cut] = "UP"
      deg[res$log2FoldChange < -log2(fc.cut) & res$padj < pval.cut] = "DN"
      res$DEG.ID = deg;
      
      # add percent expressed
      m = GetAssayData(seur[[assay.name]],slot = "counts");
      
      cc = which(colnames(m) %in% modules[gsub("^cluster.id","",cls.var[j])][[1]])
      cc.not = which(!(colnames(m) %in% modules[gsub("^cluster.id","",cls.var[j])][[1]]))
      m.in = m[match(res$gene.name,rownames(m)),cc]
      m.not = m[match(res$gene.name,rownames(m)),cc.not]
      
      res$pct.case = rowSums(m.in > 1E-320,na.rm = TRUE)/ncol(m)
      res$pct.ctrl = rowSums(m.not > 1E-320,na.rm = TRUE)/ncol(m)
      #res$parent.module = rep(vr[i],nrow(res))
      return(res)
    }
    names(cluster.res) = gsub("cluster.id","",cls.var)
    #stopCluster(cl)
    gc()
  }
  return(cluster.res)
  
}

pseudobulk_markers <- function(seu,
modules,htbl,root.module,
# seurat marker params
pval.cut = 0.05,fc.cut = 1.2,
assay.name = "SCT",group.var = "individualID",
# parallelize param
do.par = TRUE,n.cores = 8)
{
	require(Matrix.utils)
	require(DESeq2)
	require(BiocParallel)
	require(dplyr)
	require(magrittr)
	
	vr = c(root.module)
	niter = 1
	marker.res = DataFrame()
	do.run = TRUE
	while (do.run)
	{
		vr.new = c()
		cat(paste0("iter. #",niter,"\n"))
		cat(paste0("- root modules:",paste(vr,collapse = ","),"\n"))
		for (i in 1:length(vr))
		{
			mids = subset(htbl,module.parent == vr[i])$module.id
			if (length(mids) > 1)
			{
				cat(paste0("-- ",vr[i]," child modules:",paste(mids,collapse = ","),"\n"))
				mr = modules[[which(names(modules) == vr[i])]]
				mc = modules[match(mids,names(modules))]
				seur = seu[,match(mr,colnames(seu))];
				
				if (length(mc) > 1)
				{
				  ## First, aggregate cell-level to sample-level data
				  cat("-- Pseudobulk matrix generation...\n")
				  
				  pseudo.dat = create_pseudobulk_per_cluster(seur,mods = mc,group.var = group.var,assay.name = assay.name)
				  cluster.res = find_pseudobulk_markers(count.mat = pseudo.dat$pseudobulk,meta = pseudo.dat$meta,
				                                        seur = seur,modules = mc,
				                                        assay.name = assay.name,
				                                        pval.cut = pval.cut,fc.cut = fc.cut,
				                                        do.par = do.par,n.cores = n.cores)
				  # update to marker results
				  marker.res = rbind(marker.res,do.call('rbind.DataFrame',cluster.res))
				  
				}
				
				rm(seur)
				gc()
			}
			
			vr.new = c(vr.new,mids[mids %in% htbl$module.parent])
		}
		
		# update root directory
		if (length(vr.new) > 0)
		{
			do.run = TRUE
			vr = vr.new
			niter = niter + 1
			rm(vr.new)
		}else{
			do.run = FALSE
		}
	}

	return(marker.res)
}


### marker dot plot in hierarchy format
plot_hierarchy_dot <- function(seu,htbl,modules,celltype.markers,assay.name = "RNA",pct.mark = 0.1,
plot.hierarchy = TRUE,filename = "celltype_markers.hierarchy_dot.png",width.per.plot = 2500,height.per.plot = 3500)
{
	cat("Calculate marker statistics...\n")
	# get hierarchy and modules
	#htbl = clus.results$module.table[,c("cluster.name","parent.cluster")]
	#colnames(htbl) = c("module.id","module.parent")
	#modules = clus.results$modules;modules = c(list(M0 = colnames(seu)),modules)

	# get normalized data
	m = GetAssayData(seu,slot = "data",assay = assay.name)

	# per module, calculate percent of expressed genes, and average expressions
	#mids = setdiff(names(modules),"M0")
	mids = union(htbl[[1]],htbl[[2]])
	total.data = data.frame()

	for (i in 1:length(mids))
	{
		mp = subset(htbl,module.id == mids[i])$module.parent
		
		if (length(mp) > 0)
		{
		  # get parental module count
		  rr = which(rownames(seu@assays[assay.name][[1]]@counts) %in% celltype.markers)
		  cc.p = which(colnames(seu) %in% modules[mp][[1]] & !(colnames(seu) %in% modules[mids[i]][[1]]))
		  cc.i = which(colnames(seu) %in% modules[mids[i]][[1]])
		  
		  cnt.p = rowSums(seu@assays[assay.name][[1]]@counts[rr,cc.p] > 1E-320,na.rm = TRUE)
		  cnt.i = rowSums(seu@assays[assay.name][[1]]@counts[rr,cc.i] > 1E-320,na.rm = TRUE)
		  
		  # get average expressions
		  rr = which(rownames(m) %in% celltype.markers)
		  expr.p = rowMeans(m[rr,cc.p],na.rm = TRUE)
		  expr.i = rowMeans(m[rr,cc.i],na.rm = TRUE)
		  
		  # combine data
		  mrks = union(names(cnt.p),names(expr.p))
		  ct = names(celltype.markers)[match(mrks,celltype.markers)]
		  
		  dat = data.frame(module.id = rep(mids[i],length(mrks)),module.parent = rep(mp,length(mrks)),
		                   marker.gene = mrks,celltype = ct,
		                   pct.module = cnt.i[match(mrks,names(cnt.i))]/length(cc.i),
		                   pct.others = cnt.p[match(mrks,names(cnt.p))]/length(cc.p),
		                   avg.expr.module = expr.i[match(mrks,names(expr.i))],
		                   avg.expr.parent = expr.p[match(mrks,names(expr.p))])
		  
		  total.data = rbind.data.frame(total.data,dat)
		  rm(mp,cc.p,cc.i,cnt.p,cnt.i,rr,expr.p,expr.i,mrks,ct,dat)
		}
		
	}

	# add M0 stats
	rr = which(rownames(seu@assays[assay.name][[1]]@counts) %in% celltype.markers)
	cc = which(colnames(seu) %in% modules$M0)
	pct.o = rowSums(seu@assays[assay.name][[1]]@counts[rr,cc] > 1E-320,na.rm = TRUE)/length(cc)
	expr.o = rowMeans(m[rownames(m) %in% celltype.markers,colnames(seu) %in% modules$M0],na.rm = TRUE)
	mrks = union(names(pct.o),names(expr.o))

	dat.o = data.frame(module.id = rep("M0",length(mrks)), marker.gene = mrks,pct.module = pct.o[match(mrks,names(pct.o))],avg.expr.module = expr.o[match(mrks,names(expr.o))])

	### try to plot out marker expression by hierarchy format per marker gene
	plist = NULL
	if (plot.hierarchy)
	{
		cat("Plotting the hierarchy with dots...\n")
		require(ggraph)
		require(igraph)
		require(ggplot2)
		require(ggrepel)

		# create the basis graph object
		g = graph.data.frame(htbl[,c(2,1)],directed = TRUE)
		pobj = ggraph(graph = g, layout = 'dendrogram', circular = TRUE) + geom_edge_diagonal()
		
		# create plot per celltype marker
		plist = vector("list",length(celltype.markers));names(plist) = celltype.markers
		for (i in 1:length(celltype.markers))
		{
			dat = rbind.data.frame(subset(dat.o,marker.gene == celltype.markers[i]),
			subset(total.data,marker.gene == celltype.markers[i])[,c("module.id","marker.gene","pct.module","avg.expr.module")])
			
			# append expression features to node data
			node.dat = cbind.data.frame(pobj$data,dat[match(pobj$data$name,dat$module.id),-1])
			
			# find modules to mark:
			mid.lab = subset(total.data,pct.module > pct.mark & pct.module > subset(dat.o,marker.gene == celltype.markers[i])$pct.module & avg.expr.module > subset(dat.o,marker.gene == celltype.markers[i])$avg.expr.module &
			pct.module > pct.others & marker.gene == celltype.markers[i] & avg.expr.module >= avg.expr.parent)$module.id
			lab.dat = subset(node.dat,name %in% mid.lab)
			
			mrk.obj = pobj + geom_point(data = node.dat,aes(x = x,y = y,size = pct.module,colour = avg.expr.module)) + theme_bw() + 
			scale_colour_gradient2(low = "grey",mid = "white",high = "red") + 
			labs(title = paste0(names(celltype.markers)[i],":",celltype.markers[i]) )+ 
			guides(size = guide_legend(title = "%. Expressed",title.position="top", title.hjust = 0.5),colour = guide_colorbar(title = "Average\nExpresion",title.position="top", title.hjust = 0.5)) + 
			theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
				  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text = element_blank(),
				  legend.position="bottom",legend.box="horizontal")
				  
			if (nrow(lab.dat) > 0)
			{
				mrk.obj = mrk.obj + geom_text_repel(dat = lab.dat,aes(x = x,y = y,label = name))
			}
			
			plist[[i]] = mrk.obj;
			rm(node.dat,lab.dat,mrk.obj,mid.lab,dat)
	 
		}
		
		cat(paste0("output plots to:",filename,"\n"))
		require(cowplot)
		nc = ceiling(sqrt(length(plist)))
		nr = ceiling(length(plist)/nc)
		png(file = filename,res = 500,width = nc * width.per.plot,height = nr * height.per.plot)
		plot_grid(plotlist = plist,nrow = nr,ncol = nc,align = "hv")
		dev.off()
		
	}
	
	return(list(marker.stats = total.data,plotlist = plist))
}

plot_hierarchy_dot.raw <- function(m,cnts,htbl,modules,celltype.markers,pct.mark = 0.1,
                               plot.hierarchy = TRUE,filename = "celltype_markers.hierarchy_dot.png",width.per.plot = 2500,height.per.plot = 3500)
{
  cat("Calculate marker statistics...\n")
  # get hierarchy and modules
  #htbl = clus.results$module.table[,c("cluster.name","parent.cluster")]
  #colnames(htbl) = c("module.id","module.parent")
  #modules = clus.results$modules;modules = c(list(M0 = colnames(seu)),modules)
  
  # get normalized data
  #m = GetAssayData(seu,slot = "data",assay = assay.name)
  
  # per module, calculate percent of expressed genes, and average expressions
  mids = setdiff(names(modules),"M0")
  total.data = data.frame()
  
  for (i in 1:length(mids))
  {
    mp = subset(htbl,module.id == mids[i])$module.parent
    
    if (length(mp) > 0)
    {
      # get parental module count
      rr = which(rownames(cnts) %in% celltype.markers)
      cc.p = which(colnames(cnts) %in% modules[mp][[1]] & !(colnames(cnts) %in% modules[mids[i]][[1]]))
      cc.i = which(colnames(cnts) %in% modules[mids[i]][[1]])
      
      cnt.p = rowSums(cnts[rr,cc.p] > 1E-320,na.rm = TRUE)
      cnt.i = rowSums(cnts[rr,cc.i] > 1E-320,na.rm = TRUE)
      
      # get average expressions
      rr = which(rownames(m) %in% celltype.markers)
      expr.p = rowMeans(m[rr,cc.p],na.rm = TRUE)
      expr.i = rowMeans(m[rr,cc.i],na.rm = TRUE)
      
      # combine data
      mrks = union(names(cnt.p),names(expr.p))
      ct = names(celltype.markers)[match(mrks,celltype.markers)]
      
      dat = data.frame(module.id = rep(mids[i],length(mrks)),module.parent = rep(mp,length(mrks)),
                       marker.gene = mrks,celltype = ct,
                       pct.module = cnt.i[match(mrks,names(cnt.i))]/length(cc.i),
                       pct.others = cnt.p[match(mrks,names(cnt.p))]/length(cc.p),
                       avg.expr.module = expr.i[match(mrks,names(expr.i))],
                       avg.expr.parent = expr.p[match(mrks,names(expr.p))])
      
      total.data = rbind.data.frame(total.data,dat)
      rm(mp,cc.p,cc.i,cnt.p,cnt.i,rr,expr.p,expr.i,mrks,ct,dat)
    }
    
  }
  
  # add M0 stats
  rr = which(rownames(cnts) %in% celltype.markers)
  cc = which(colnames(cnts) %in% modules$M0)
  pct.o = rowSums(cnts[rr,cc] > 1E-320,na.rm = TRUE)/length(cc)
  expr.o = rowMeans(m[rownames(m) %in% celltype.markers,colnames(cnts) %in% modules$M0],na.rm = TRUE)
  mrks = union(names(pct.o),names(expr.o))
  
  dat.o = data.frame(module.id = rep("M0",length(mrks)), marker.gene = mrks,pct.module = pct.o[match(mrks,names(pct.o))],avg.expr.module = expr.o[match(mrks,names(expr.o))])
  
  ### try to plot out marker expression by hierarchy format per marker gene
  plist = NULL
  if (plot.hierarchy)
  {
    cat("Plotting the hierarchy with dots...\n")
    require(ggraph)
    require(igraph)
    require(ggplot2)
    require(ggrepel)
    
    # create the basis graph object
    g = graph.data.frame(htbl[,c(2,1)],directed = TRUE)
    pobj = ggraph(graph = g, layout = 'dendrogram', circular = TRUE) + geom_edge_diagonal()
    
    # create plot per celltype marker
    plist = vector("list",length(celltype.markers));names(plist) = celltype.markers
    for (i in 1:length(celltype.markers))
    {
      dat = rbind.data.frame(subset(dat.o,marker.gene == celltype.markers[i]),
                             subset(total.data,marker.gene == celltype.markers[i])[,c("module.id","marker.gene","pct.module","avg.expr.module")])
      
      # append expression features to node data
      node.dat = cbind.data.frame(pobj$data,dat[match(pobj$data$name,dat$module.id),-1])
      
      # find modules to mark:
      mid.lab = subset(total.data,pct.module > pct.mark & pct.module > subset(dat.o,marker.gene == celltype.markers[i])$pct.module & avg.expr.module > subset(dat.o,marker.gene == celltype.markers[i])$avg.expr.module &
                         pct.module > pct.others & marker.gene == celltype.markers[i] & avg.expr.module >= avg.expr.parent)$module.id
      lab.dat = subset(node.dat,name %in% mid.lab)
      
      mrk.obj = pobj + geom_point(data = node.dat,aes(x = x,y = y,size = pct.module,colour = avg.expr.module)) + theme_bw() + 
        scale_colour_gradient2(low = "grey",mid = "white",high = "red") + 
        labs(title = paste0(names(celltype.markers)[i],":",celltype.markers[i]) )+ 
        guides(size = guide_legend(title = "%. Expressed",title.position="top", title.hjust = 0.5),colour = guide_colorbar(title = "Average\nExpresion",title.position="top", title.hjust = 0.5)) + 
        theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
              panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text = element_blank(),
              legend.position="bottom",legend.box="horizontal")
      
      if (nrow(lab.dat) > 0)
      {
        mrk.obj = mrk.obj + geom_text_repel(dat = lab.dat,aes(x = x,y = y,label = name))
      }
      
      plist[[i]] = mrk.obj;
      rm(node.dat,lab.dat,mrk.obj,mid.lab,dat)
      
    }
    
    if (!is.null(filename))
    {
      cat(paste0("output plots to:",filename,"\n"))
      require(cowplot)
      nc = ceiling(sqrt(length(plist)))
      nr = ceiling(length(plist)/nc)
      png(file = filename,res = 500,width = nc * width.per.plot,height = nr * height.per.plot)
      plot_grid(plotlist = plist,nrow = nr,ncol = nc,align = "hv")
      dev.off()
      
    }
    
  }
  
  return(list(marker.stats = total.data,plotlist = plist))
}

# pie chart to show sample makeup
overlap_modules_wt_categ <- function(mods,vec)
{
  tpe = setdiff(unique(vec),NA)
  mat = do.call('rbind',lapply(mods,function(x,y,t) {out = table(y[names(y) %in% x])/length(x);out[match(tpe,names(out))]},y = vec,t = tpe))
  mat[is.na(mat)] = 0;
  #rownames(mat) = gsub("^c1_","M",rownames(mat));
  colnames(mat) = tpe
  
  if (any(rowSums(mat,na.rm = TRUE) < 1))
  {
    mat = cbind(mat,1-rowSums(mat,na.rm = TRUE))
    colnames(mat)[ncol(mat)] = "unassigned"
  }
  return(mat)
}


plot_hierarchy_pie <- function(htbl,modules,vec,edge_colour = "black",colmap = NULL)
{
  # create the basis graph object
  g = graph.data.frame(htbl[,c(2,1)],directed = TRUE)
  pobj = ggraph(graph = g, layout = 'dendrogram', circular = TRUE) + geom_edge_diagonal(edge_colour = edge_colour)
  
  # make piechart
  library(scatterpie)
  categ.mat = overlap_modules_wt_categ(mods = modules,vec = vec )
  node.dat = cbind.data.frame(pobj$data,as.data.frame(categ.mat[match(pobj$data$name,rownames(categ.mat)),]))
  
  if (is.null(colmap))
  {
    library(RColorBrewer)
    colmap = colorRampPalette(brewer.pal(11,"Spectral"))(ncol(categ.mat))
    names(colmap) = colnames(categ.mat)
    if (any(names(colmap) == "unassigned")) colmap["unassigned"] = "grey"
  }
  
  pobj = pobj + geom_scatterpie(aes_string(x="x", y = "y"), data=node.dat,
                                cols=colnames(categ.mat), color=NA) + scale_fill_manual(values = colmap) + coord_equal() +
    theme_bw() + 
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.title = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text = element_blank(),
          legend.position="bottom",legend.box="horizontal")
  return(pobj)
}


### reduced dimension plot with customized labeling and category to color.
plot_dimred_manual <- function(seu.obj,dimred = "umap",category.name,cls.col = NULL,add.labels = TRUE)
{
	require(RColorBrewer)
	require(grid)
	require(ggplot2)
	
	xy = Embeddings(seu.obj,dimred)
	nc = nrow(xy)
	cell.names = rownames(xy)
	cls.col = NULL
	
	cell.feat = seu.obj@meta.data[[which(colnames(seu.obj@meta.data) == category.name)]]

	pdata = data.frame(x = xy[,1],y = xy[,2],cell.feature = cell.feat)

	mods = split(colnames(seu.obj),factor(cell.feat))

	if (is.null(cls.col))
	{
		if (length(mods) <= 11)
		{
			cls.col = brewer.pal(length(mods),"Spectral");
			names(cls.col) = names(mods)
		}else{
			getPalette = colorRampPalette(brewer.pal(9, "Set1"))
			cls.col = getPalette(length(mods))
			names(cls.col) = names(mods)
		}
	}

	pobj = ggplot() + geom_point(data = pdata,aes(colour = cell.feature,x = x,y = y)) + 
	scale_colour_manual(values = cls.col) + theme_bw() 
		
	# add cluster label
	if (add.labels)
	{
		cls.coord = get_module_label_coord(coord = xy,mods = mods)
		pobj = pobj + geom_text(data = data.frame(cluster = rownames(cls.coord),x = cls.coord[,1],y = cls.coord[,2]),aes(x = x,y = y,label = cluster))
	}
	return(pobj)
}

plot_dimred_manual.raw <- function(xy,cell.feat,alpha.val = 0.1,cls.col = NULL,add.labels = TRUE,label.size = 5,is.repel = FALSE,is.root = NULL)
{
  require(RColorBrewer)
  require(grid)
  require(ggplot2)
  
  #xy = Embeddings(seu.obj,dimred)
  nc = nrow(xy)
  cell.names = rownames(xy)
  cls.col = NULL
  
  if (is.null(is.root)) is.root = rep(TRUE,nrow(xy))
  pdata = data.frame(x = xy[,1],y = xy[,2],cell.feature = cell.feat,is.root = is.root)
  
  mods = split(names(cell.feat),factor(cell.feat))
  
  if (is.null(cls.col))
  {
    if (length(mods) <= 11)
    {
      cls.col = brewer.pal(length(mods),"Spectral");
      names(cls.col) = names(mods)
    }else{
      getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      cls.col = getPalette(length(mods))
      names(cls.col) = names(mods)
    }
  }
  
  pobj = ggplot() + geom_point(data = pdata,aes(colour = cell.feature,x = x,y = y),alpha = alpha.val) + 
    scale_colour_manual(values = cls.col) + 
    theme_bw() 
  
  # add cluster label
  if (add.labels)
  {
    cls.coord = get_module_label_coord(coord = xy,mods = mods)
    if (!is.repel) pobj = pobj + geom_text(data = data.frame(cluster = rownames(cls.coord),x = cls.coord[,1],y = cls.coord[,2]),aes(x = x,y = y,label = cluster),size = label.size)
    if (is.repel) 
    {
    	require(ggrepel)
    	pobj = pobj + geom_text_repel(data = data.frame(cluster = rownames(cls.coord),x = cls.coord[,1],y = cls.coord[,2]),aes(x = x,y = y,label = cluster),size = label.size)
    }
  }
  return(pobj)
}

get_rooted_membership <- function(mtbl,modules,pid,bg)
{
  vec = rep(NA,length(bg))
  names(vec) = bg
  
  mods = modules[subset(mtbl,parent.cluster == pid)$cluster.name]
  for (i in 1:length(mods)) vec[bg %in% mods[[i]]] = names(mods)[i]
  return(vec)
}
