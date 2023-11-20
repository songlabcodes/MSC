# get coordinates for centered labels 
get_module_label_coord <- function(coord,mods)
{
	do.call('rbind',lapply(mods,function(x,y) apply(y[match(x,rownames(y)),],2,function(q) median(q,na.rm = T)),y = coord))
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

#' @title Hierarchical representation of marker dot plots  
#' @name plot_hierarchy_dot.raw
#' @description Create circular hierarchy tree to visualize dot plot for desired marker genes, and calculate the marker statistics for provided genes.
#' @param m log-normalized gene expression matrix to calculate the average expression levels.  
#' @param cnts count matrix to calculate the %. of cells expressing genes. 
#' @param htbl A two column data.frame table. The first column is cluster name, and the second column is its parent cluster name.  
#' @param modules A named list object. Each is the character vector of cells from a cell cluster. 
#' @param celltype.markers A named character vector for the list of genes to create the dotplot. The names specify cell type names. 
#' @param plot.hierarchy A logical. If FALSE, only marker statistics will be returned. 
#' @param filename A character. If provided, a .png file will be outputed with the filename. 
#' @param width.per.plot A numeric. The width per plot.
#' @param height.per.plot A numeric. The height per plot.  
#' @return a list of data.frame for the marker statistics, and list of ggplot objects holding the plots.
#' @examples 
#' # See Vignettes. 
#' @export
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
      
      cnt.p = Matrix::rowSums(cnts[rr,cc.p] > 1E-320,na.rm = TRUE)
      cnt.i = Matrix::rowSums(cnts[rr,cc.i] > 1E-320,na.rm = TRUE)
      
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
  pct.o = Matrix::rowSums(cnts[rr,cc] > 1E-320,na.rm = TRUE)/length(cc)
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

#' @title Hierarchical representation of modules with piecharts  
#' @name plot_hierarchy_pie
#' @description Create circular hierarchy tree to visualize the cluster hierarchy, with pie chart to show cell cluster composition 
#' @param htbl A two column data.frame table. The first column is cluster name, and the second column is its parent cluster name.  
#' @param modules A named list object. Each is the character vector of cells from a cell cluster. 
#' @param vec A factor vector, containing different groups of cells. 
#' @param edge_colour edge colour of pie chart.
#' @param colmap a named vector to assign colors to each group in vec. 
#' @return a ggplot object
#' @examples 
#' data(pbmc_8k_msc_results)
#' htbl = pbmc_8k_msc_results$pruned.table[,c("cluster.name","parent.cluster")]
#' modules = pbmc_8k_msc_results$modules
#' bg = Reduce("union",modules)
#' vec = sample(LETTERS[1:10],length(bg),replace = TRUE);names(vec) = bg;vec = factor(vec)
#' obj = plot_hierarchy_pie(htbl,modules,vec)
#' @export
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

#' @title Generate dimension reduction plot.   
#' @name plot_dimred_manual.raw
#' @description  Generated ggplot object to show 2-dimensional dimension reduction plot 
#' @param xy Two column matrix, where each row is the cartesian coordinate for a cell. 
#' @param cell.feat A named factor vector for cell group assignments. 
#' @param alpha.val A numeric from 0 to 1 to specific dot transparency. 
#' @param cls.col A named character vector to specify colors for the groups in "cell.feat".
#' @param add.labels A logical. If TRUE, adds labels to the groups. 
#' @param label.size If add.labels = TRUE, then label.size specifies the label font size in the plot. 
#' @param is.repel A logical. If TRUE, uses [geom_text_repel()] from ggrepel package to add the labels.
#' @param is.root A logical vector. Each specify if each cell belongs to the parent cluster of interest. 
#' @return Returns a named vector of cluster membership.
#' @examples 
#' cell.names = letters[1:10]
#' xy =  matrix(runif(20,0,1),ncol = 2);rownames(xy) = cell.names
#' cell.feat= factor(c(rep("A",5),rep("B",5)));names(cell.feat) = cell.names
#' plot_dimred_manual.raw(xy,cell.feat,alpha.val = 1)
#' @export

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

#' @title Get child cluster membership vector  
#' @name get_rooted_membership
#' @description  Anchored onto a parent cluster, child cluster membership vector is calculated 
#' @param mtbl A data.frame object, for the cluster information table. Can be obtained from module.table or pruned.table slot in output from [iterative_clustering.par()]. 
#' @param modules A named list object. Each is the character vector of cells from a cell cluster. 
#' @param pid A character value for the parent cluster. The respective child cluster will be identified from mtbl. 
#' @param bg A character vector containing the full list of cell names.
#' @return Returns a named vector of cluster membership.
#' @examples 
#' data(pbmc_8k_msc_results)
#' library(igraph)
#' mem.vec = get_rooted_membership(mtbl = pbmc_8k_msc_results$pruned.table,modules = pbmc_8k_msc_results$modules,pid = "M2",bg = V(pbmc_8k_msc_results$cell.network)$name)
#' table(mem.vec)
#' @export

get_rooted_membership <- function(mtbl,modules,pid,bg)
{
  vec = rep(NA,length(bg))
  names(vec) = bg
  
  mods = modules[subset(mtbl,parent.cluster == pid)$cluster.name]
  for (i in 1:length(mods)) vec[bg %in% mods[[i]]] = names(mods)[i]
  return(vec)
}
