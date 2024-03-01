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

ij_to_lin <- function(rr,cc,n) n*(cc-1) + rr

convert_mat_to_edgelist <- function(m)
{
  ## quickly convert matrix to three column list
  ij = do.call('rbind',lapply(1:(ncol(m)-1),function(x,n) {rr = rep(x);cc = (x+1):n;cbind(rr,cc)},n = ncol(m)));
  colnames(ij) = c("ri","ci")
  w = m[ij_to_lin(rr = ij[,1],cc = ij[,2],n = ncol(m))]
  cbind(ij,w)
}

get_jaccard <- function(count.dat,min.count = 1)
{
  mat.dat.bin = count.dat
  mat.dat.bin[count.dat >= min.count] = 1
  mat.dat.bin[count.dat < min.count] = 0
  
  nzero.cell = colSums(mat.dat.bin)
  
  # get jaccard index
  ints = as.matrix(t(mat.dat.bin) %*% mat.dat.bin)
  norm.factor = matrix(rep(nzero.cell,ncol(mat.dat.bin)),nrow = ncol(mat.dat.bin))
  jac = ints/(norm.factor + t(norm.factor) - ints)
  dimnames(jac) = dimnames(ints)
  return(jac)
}

adj_to_el <- function(adj,is.upper.tri = TRUE)
{
  if (is.upper.tri) 
  {
    ij = which(adj > 1E-320 & upper.tri(adj),arr.ind = TRUE)
  }else{
    ij = which(adj > 1E-320,arr.ind = TRUE)
  }
  return(ij)
}

get_local_embedding_per_node <- function(v,mat,min.k = 3,max.k = 30,n_skip = 3,
                                         dist.func = function(a, b) sqrt(sum((a - b)^2)),
                                         is.decreasing = FALSE,verbose = FALSE)
{
  convert_mat_to_edgelist <- function(m)
  {
    ## quickly convert matrix to three column list
    ij = do.call('rbind',lapply(1:(ncol(m)-1),function(x,n) {rr = rep(x);cc = (x+1):n;cbind(rr,cc)},n = ncol(m)));
    colnames(ij) = c("ri","ci")
    w = m[ij_to_lin(rr = ij[,1],cc = ij[,2],n = ncol(m))]
    cbind(ij,w)
  }
  ij_to_lin <- function(rr,cc,n) n*(cc-1) + rr
	idx = setdiff(1:ncol(mat),v)
	
	#rval = apply(mat[,idx],2,function(x,y) cor(x,y,method = method),y = mat[,v])
	#rval = cor(x = mat[,v],y = mat[,idx],method = method,use = "na.or.complete");
	rval = apply(mat[,idx],2,function(x,y) dist.func(x,y),y = mat[,v])
	if (!is.vector(rval)) rval = as.vector(rval)
	
	if (max.k > length(rval)) max.k = length(rval)
	
	# get indices for max.k 
	ii = order(rval,decreasing = is.decreasing)
	rval = rval[ii]
	idx =idx[ii]
	max.k.idx = idx[1:max.k]
	max.k.rval = rval[1:max.k]

	# Now, calculate correlation matrix within max.k 
	top.idx = c(v,max.k.idx)
	top.rho = matrix(0,nrow = length(top.idx),ncol = length(top.idx));
	rownames(top.rho) = colnames(top.rho) = colnames(mat)[top.idx]
	
	for (i in 1:(length(top.idx)-1))
	{
	  for (j in (i + 1):length(top.idx))
	  {
	    val = dist.func(mat[,top.idx[i]],mat[,top.idx[j]])
	    top.rho[i,j] = top.rho[j,i] = val 
	  }
	}
	
	##### initiate while loop parameters
	cur_link = 1
	n_o = min.k;
	skipped = 0
	do.run = TRUE
	nm_o = colnames(top.rho)[1]
	knn.f = 0
	## repeat while loop until saturation
	if (verbose) cat("knn=")

	while (do.run)
	{
		  if (verbose) cat(paste((n_o),",",sep = ""))
		  
		  # create edgelist for top n_o nodes
		  ord.idx = 1:(n_o + 1);
		  ijw_n = convert_mat_to_edgelist(m = top.rho[ord.idx,ord.idx])
		  ijw_n <- ijw_n[order(ijw_n[,3],decreasing = is.decreasing),]
		  ijw_n = data.frame(row = colnames(top.rho)[ord.idx[ijw_n[,1]]],col = colnames(top.rho)[ord.idx[ijw_n[,2]]],weight = ijw_n[,3])
		  
		  sink("tmp.mst")
		  el <- calculate.PFN(ijw_n,doPar = FALSE)
		  sink(NULL)
		  el <- subset(el,row == nm_o)
		  
		  if (nrow(el) > cur_link)
		  {
  			ijw_f <- el;
  			cur_link <- nrow(el);
		  }else{
			  skipped = skipped + 1
		  }
		  #ij_o <- ij_n
		  #ijw_o <- ijw_n
		  
		  n_o = n_o +1 
		  if (skipped >= n_skip) 
		  {
			  do.run = FALSE
			  if (verbose) cat("\n")
			  knn.f = n_o-1
		  }
		  
		  if (n_o > max.k)
		  {
  			do.run = FALSE
  			if (verbose) cat("\n")
  			knn.f = n_o-1
		  }
		  rm(ijw_n,el)
	}

	return(ijw_f)
}


get_local_embedding_par <- function(mat,n.cores = 4,n_skip = 3,min.k = 4,max.k = 30,
                                    dist.func = function(a, b) sqrt(sum((a - b)^2)),is.decreasing = FALSE,
                                    verbose = FALSE)
{
  ### Set up parallel backend
  require(doParallel)
  require(foreach)
  
  if (getDoParWorkers() == 1 & n.cores > 1)
  {
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    # check how many workers are there
    cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
  }
  
  # get dim
  nv = ncol(mat)
  
  # split jobs 
  dn <- ceiling(nv/n.cores)
  fact <- factor(do.call('c',lapply(1:n.cores,function(i,dn) rep(i,dn),dn = dn))[1:nv])
  idx.lst = split(1:nv,fact)
  
  cat("Job split into:\n");
  print(sapply(idx.lst,length))
  
  # get local embedding per node
  res <- foreach(ii = idx.lst,.packages = "MEGENA",.export = c("get_local_embedding_per_node","convert_mat_to_edgelist","ij_to_lin")) %dopar% {
    
    res.ijw = data.frame()
    knn.f = c()
    
    for (i in ii)
    {
      ijw.v = get_local_embedding_per_node(v = i,mat,min.k = min.k,max.k = max.k,n_skip = n_skip,
                                           dist.func = dist.func,is.decreasing = is.decreasing,verbose = verbose)
      res.ijw = rbind.data.frame(res.ijw,ijw.v)
      rm(out)
    }
    
    return(res.ijw)
  }
  
  out = do.call('rbind.data.frame',res)
  out = subset(out,as.character(row) != as.character(col))
  return(out)
}

# this version works on sparse count matrix to speed up the process
calculate_LEN.v3 <- function(mat,bcnt,
                             do.par = TRUE,n.cores = 4,
                             dist.func = function(a, b) cor(x = a,y = b,method = "pearson"),
                             is.decreasing = TRUE,min.k = 3,max.k = 30,n_skip = 3)
{
  ##### Explore if we can further speed up LENA construction by subsetting cells that share 
  
  
  ## Set up parallel backend
  require(doParallel)
  require(foreach)
  
  if (getDoParWorkers() == 1 & n.cores > 1 & do.par)
  {
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    # check how many workers are there
    cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
  }
  
  
  ## subset local search size that substantially share expressed genes
  #cat(paste0("Subsetting cells sharing similar list of expressed genes...\n"))
  ng.per.cell = colSums(bcnt,na.rm = T)
  
  # screen out cells with no gene expressions
  if (any(ng.per.cell < 1E-320))
  {
    ii = ng.per.cell > 1E-320
    m = m[,ii]
    bcnt = bcnt[,ii]
    ng.per.cell = ng.per.cell[ii]
    rm(ii)
  }
  gc()
  
  ## split jobs
  dn = ceiling(ncol(mat)/n.cores)
  
  cat(paste0("- #. jobs per node:",dn,"\n"))
  fact = factor(do.call('c',lapply(1:n.cores,function(i) rep(i,dn)))[1:ncol(mat)])
  idx.lst = split(1:ncol(mat),fact)
  
  # Start the clock!
  ptm <- proc.time()
  
  res <- foreach(ii = idx.lst,.packages = "MEGENA",
                 .export = c("get_local_embedding_per_node","convert_mat_to_edgelist","ij_to_lin","ng.per.cell","min.k","max.k")) %dopar% {
                   res.ijw = data.frame()
                   knn.f = c()
                   
                   for (n in ii)
                   {
                     # calculate pairs to run the analysis by checking shared gene expressions
                     vec = as.vector(as.vector(bcnt[,n]) %*% bcnt)
                     jr = vec/(ng.per.cell - vec + ng.per.cell[n])
                     
                     qc = quantile(jr,probs = (length(jr) - 2*max.k)/length(jr))
                     jj = which(jr >= qc)
                     
                     if (length(jj) > min.k)
                     {
                       ijw.v = get_local_embedding_per_node(v = 1,mat = mat[,c(n,setdiff(jj,n))],min.k = min.k,max.k = min(c(length(jj),max.k)),n_skip = n_skip,
                                                            dist.func = dist.func,is.decreasing = is.decreasing,verbose = FALSE)
                       res.ijw = rbind.data.frame(res.ijw,ijw.v)
                     }
                   }
                   return(res.ijw)
                 }
  
  # Stop the clock
  cat("- runtime to calculate LEN per node:\n")
  print(proc.time() - ptm)
  
  el = do.call('rbind.data.frame',res)
  el = subset(el,as.character(row) != as.character(col) & weight > 1E-320) # ensure positivity in weights 
  
  # symmetrize adjacency
  adj = sparseMatrix(i = match(el$row,colnames(mat)),j = match(el$col,colnames(mat)),x = el$weight,dims = c(ncol(mat),ncol(mat)),dimnames = list(colnames(mat),colnames(mat)))
  adj = adj + t(adj)		
  
  ### get fractional ranking
  ptm <- proc.time()
  cat("- Check duplicates...\n")
  is.mutual = is.duplicate = rep(FALSE,nrow(el))
  for (i in 1:nrow(el))
  {
    if ((i %% 1000) == 0) cat(paste0("processed:",i,"/",nrow(el),"=",signif(i/nrow(el),3)*100,"%\n"))
    ii = c()
    if (!is.duplicate[i]) ii = which(el$row == el$col[i] & el$col == el$row[i])
    if (length(ii) > 0) 
    {
      is.mutual[c(i,ii)] = TRUE
      is.duplicate[ii] = TRUE
    }
  }
  
  # get jaccard index
  cat("- Calculate mutual neighbor ratios...\n")
  jac <-rep(NA,nrow(el))
  for (i in 1:nrow(el))
  {
    if ((i %% 1000) == 0) cat(paste0("processed:",i,"/",nrow(el),"=",signif(i/nrow(el),3)*100,"%\n"))
    
    vi = which(rownames(adj) == el$row[i])
    vj = which(rownames(adj) == el$col[i])
    
    veci = adj[vi,] > 1E-320
    vecj = adj[vj,] > 1E-320
    veci[c(vi,vj)] = TRUE
    vecj[c(vi,vj)] = TRUE
    jac[i] = sum(veci & vecj)/sum(veci | vecj,na.rm = TRUE)
  }			  
  el$Mutual.Neighbor.Ratio = jac;
  
  el$is.mutual = is.mutual
  el$is.duplicate = is.duplicate
  
  cat("- runtime to calculate mutual neighbor ratios:\n")
  print(proc.time() - ptm)
  
  return(el)
}
# this version works on normalized expression only
calculate_LEN.v4 <- function(mat,
                             do.par = TRUE,n.cores = 4,
                             dist.func = function(a, b) cor(x = a,y = b,method = "pearson"),
                             is.decreasing = TRUE,min.k = 3,max.k = 30,n_skip = 3)
{
  ## Set up parallel backend
  require(doParallel)
  require(foreach)
  
  if (getDoParWorkers() == 1 & n.cores > 1 & do.par)
  {
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    # check how many workers are there
    cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
  }
  
  ## split jobs
  dn = ceiling(ncol(mat)/n.cores)
  
  cat(paste0("- #. jobs per node:",dn,"\n"))
  fact = factor(do.call('c',lapply(1:n.cores,function(i) rep(i,dn)))[1:ncol(mat)])
  idx.lst = split(1:ncol(mat),fact)
  
  # Start the clock!
  ptm <- proc.time()
  
  res <- foreach(ii = idx.lst,.packages = "MEGENA",
                 .export = c("get_local_embedding_per_node","convert_mat_to_edgelist","ij_to_lin","min.k","max.k")) %dopar% {
                   res.ijw = data.frame()
  
                   knn.f = c()
                   
                   for (n in ii)
                   {
                     out = apply(mat,2,function(x,y) dist.func(x,y),y = mat[,n])
                     jj = order(out,decreasing = is.decreasing)[1:(max.k + 1)]
                     if (length(jj) > min.k)
                     {
                       ijw.v = get_local_embedding_per_node(v = 1,mat = mat[,c(n,setdiff(jj,n))],min.k = min.k,max.k = min(c(length(jj),max.k)),n_skip = n_skip,
                                                            dist.func = dist.func,is.decreasing = is.decreasing,verbose = FALSE)
                       res.ijw = rbind.data.frame(res.ijw,ijw.v)
                     }
                   }
                   return(res.ijw)
                 }
  
  # Stop the clock
  cat("- runtime:\n")
  print(proc.time() - ptm)
  
  el = do.call('rbind.data.frame',res)
  el = subset(el,as.character(row) != as.character(col) & weight > 1E-320) # ensure positivity in weights 
  
  # symmetrize adjacency
  adj = sparseMatrix(i = match(el$row,colnames(mat)),j = match(el$col,colnames(mat)),x = el$weight,dims = c(ncol(mat),ncol(mat)),dimnames = list(colnames(mat),colnames(mat)))
  adj = adj + Matrix::t(adj)		
  
  ### get fractional ranking
  cat("- Check duplicates...\n")
  is.mutual = is.duplicate = rep(FALSE,nrow(el))
  for (i in 1:nrow(el))
  {
    if ((i %% 1000) == 0) cat(paste0("processed:",i,"/",nrow(el),"=",signif(i/nrow(el),3)*100,"%\n"))
    ii = c()
    if (!is.duplicate[i]) ii = which(el$row == el$col[i] & el$col == el$row[i])
    if (length(ii) > 0) 
    {
      is.mutual[c(i,ii)] = TRUE
      is.duplicate[ii] = TRUE
    }
  }
  
  # get jaccard index
  cat("- Calculate mutual neighbor ratios...\n")
  jac <-rep(NA,nrow(el))
  for (i in 1:nrow(el))
  {
    if ((i %% 1000) == 0) cat(paste0("processed:",i,"/",nrow(el),"=",signif(i/nrow(el),3)*100,"%\n"))
    
    vi = which(rownames(adj) == el$row[i])
    vj = which(rownames(adj) == el$col[i])
    
    veci = adj[vi,] > 1E-320
    vecj = adj[vj,] > 1E-320
    veci[c(vi,vj)] = TRUE
    vecj[c(vi,vj)] = TRUE
    jac[i] = sum(veci & vecj)/sum(veci | vecj,na.rm = TRUE)
  }			  
  el$Mutual.Neighbor.Ratio = jac;
  
  el$is.mutual = is.mutual
  el$is.duplicate = is.duplicate
  
  return(el)
}

calculate_mutual_neighbor_ratio <- function(g,vi,vj)
{
  neigh = lapply(ego(graph = g,order = 1,nodes = c(vi,vj),mode = "all"),function(x) names(x))
  nr = length(intersect(neigh[[1]],neigh[[2]]))/sum(sapply(neigh,length))
  return(nr)
}

prune_edges <- function(dat.el,quantile.cut = 0.3)
{
  go = g = graph.data.frame(dat.el,directed = FALSE)
  thresh = quantile(subset(dat.el,!is.mutual)$Mutual.Neighbor.Ratio,probs = quantile.cut)
  
  mut.el.o = subset(dat.el[order(dat.el$Mutual.Neighbor.Ratio),],Mutual.Neighbor.Ratio < thresh & !is.mutual)
  mut.el.obj = paste(mut.el.o[[1]],mut.el.o[[2]],sep = "|")
  do.remove = rep(FALSE,nrow(mut.el.o))
  for (r in 1:nrow(mut.el.o))
  {
    if ((r %% 1000) == 0) cat(paste("Pruning:",signif(r/nrow(mut.el.o)*100,3),"%\n"))
    mut.o = calculate_mutual_neighbor_ratio(g = go,vi = mut.el.o[[1]][r],vj = mut.el.o[[2]][r])
    gf = go - edges(mut.el.obj[r])
    mut.f = calculate_mutual_neighbor_ratio(g = gf,vi = mut.el.o[[1]][r],vj = mut.el.o[[2]][r])
    
    sout = shortest_paths(graph = gf,from = mut.el.o[[1]][r],to = mut.el.o[[2]][r])
    if (length(sout$vpath[[1]]) > 0 & mut.o > 0 & mut.f < mut.o) 
    {
      do.remove[r] = TRUE
      go = gf;
      rm(gf)
    }
    
  }
  
  return(go)
}

#' @title Perform locally embedded network (LEN) to calculate cell-cell similarity network
#' @name generate_cell_network.wt.loess
#' @description  Perform locally embedded network (LEN) to calculate cell-cell similarity network. 
#' @param mat A log-normalized single-cell expression matrix (rows: genes, columns: cell).
#' @param bcnt A binary matrix of genes (row) by cells (column). If expressed, \code{bcnt[i,j] = 1}, otherwise 0.   
#' @param dist.func A function to calculate cell-cell similarity/dissimilarity. 
#' @param is.decreasing A logical. TRUE if dist.func is similarity function to show greater value corresponds to greater similarity. FALSE if dist.func is dissimilarity.  
#' @param n.cores An integer for the number of cores for parallelization 
#' @return 
#' A list containing the MSC results
#' \describe{
#'   \item{modules}{A list of cell clusters}
#'   \item{module.table}{A data.frame containing cell cluster statistics (sizes, compactness, intra-cluster connectivity etc).}
#'   \item{cell.network}{An igraph object containing the cell network}
#'   \item{alpha.value}{The alpha parameter value used to calculate the cluster compactness}
#'   \item{pruned.table}{Refined list of cell clusters after applying compactness and/or intra-cluster connectivity filters.}
#' }
#' @examples 
#' \donttest{
#' data(simMix1)
#' 
#' library(scater)
#' library(scran)
#' library(doParallel)
#' qc.out = process_data(simMix1,do.impute = TRUE)
#' 
#' # add reduced dim
#' reducedDim(qc.out$sce,"PCA.imputed") = calculatePCA(x = assay(qc.out$sce,"ALRA.logcounts"),subset_row = rownames(subset(qc.out$gene.data.alra,p.value < 0.05)))
#' reducedDim(qc.out$sce,"UMAP.imputed") = calculateUMAP(x = qc.out$sce,dimred = "PCA.imputed")
#' reducedDim(qc.out$sce,"tSNE.imputed") = calculateTSNE(x = qc.out$sce,dimred = "PCA.imputed")
#' 
#' plotReducedDim(qc.out$sce,dimred = "tSNE.imputed",colour_by = "phenoid")

#' # for imputed data for correlation computation 
#' dat = subset(qc.out$gene.data.alra,bio > 0 & p.value < 0.05)
#' bcnt.alra = counts(qc.out$sce)[rownames(dat),]
#' m.alra = as.matrix(assay(qc.out$sce,"ALRA.logcounts")[rownames(dat),])
#' bcnt.alra[bcnt.alra > 1E-320] = 1
#' g.len = generate_cell_network.wt.loess(mat = m.alra,bcnt = bcnt.alra,dist.func = function(a, b) cor(a,b,method ="pearson"),is.decreasing = TRUE,n.cores = 4)
#' }
#' @export
generate_cell_network.wt.loess <- function(mat,bcnt,dist.func,is.decreasing,n.cores = 8)
{
  require(rpart)
  # get network
  if (getDoParWorkers() == 1 & n.cores > 1 & !(n.cores == getDoParWorkers()))
  {
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
  }
  
  # use binarized expression matrix if present (this speeds up substantially)
  if (!is.null(bcnt))
  {
    el = calculate_LEN.v3(mat,bcnt,do.par = TRUE,n.cores = n.cores,
                          dist.func = dist.func,
                          is.decreasing = is.decreasing,min.k = 3,max.k = 30,n_skip = 3)
    
  }else{
    el = calculate_LEN.v4(mat,do.par = TRUE,n.cores = n.cores,
                          dist.func = dist.func,
                          is.decreasing = is.decreasing,min.k = 3,max.k = 30,n_skip = 3)
    
  }
  #el = calculate_LEN(mat,bcnt,n.cores = n.cores,
  #                      dist.func = dist.func,
  #                      is.decreasing = is.decreasing,min.k = 3,max.k = 30,n_skip = 3)
  
  gc()
  
  dat.el = subset(el,!is.duplicate & weight > 1E-320)
  dat.el = dat.el[order(dat.el$Mutual.Neighbor.Ratio,decreasing = TRUE),]
  
  fit.loess = loess(formula= weight ~ Mutual.Neighbor.Ratio,data = dat.el,span = 0.2)
  smoothed <- predict(fit.loess,se= FALSE)
  resid.val = residuals(fit.loess)
  sigma = sd(resid.val)
  plot(dat.el[,c("Mutual.Neighbor.Ratio","weight")])
  lines(smoothed, x=dat.el$Mutual.Neighbor.Ratio, col="red")
  lines(x=dat.el$Mutual.Neighbor.Ratio,smoothed - 2*sigma, lty=2,col = "red")
  lines(x=dat.el$Mutual.Neighbor.Ratio,smoothed + 2*sigma, lty=2,col = "red")
  
  # find which links are 
  if (is.decreasing) dat.el$is.low.weight = residuals(fit.loess) < -2*sigma
  if (!is.decreasing) dat.el$is.low.weight = residuals(fit.loess) > 2*sigma
  print(table(dat.el$is.low.weight))
  
  gf = prune_edges(subset(dat.el,!is.low.weight),quantile.cut = 0.25)
  
  return(gf)
}

###### Network visualization toolbox
find_elbow_dist <- function(x,y,cushion = 0.05,verbose = FALSE)
{
    #### Assume monotonically decreasing curve
    # get max and min
    ymx = max(y,na.rm = TRUE)
    ymn = min(y,na.rm = TRUE)
    
    # discard minimum plateau
    #ii = min(which(y <= (ymn + 1E-320)))
    #xi = x[1:ii];yi = y[1:ii]
    xi = x;yi = y
    
    ### get distance 
    # get angle with reference line
    veco = c(xi[length(xi)] - xi[1],yi[length(yi)] - yi[1])
    veci = cbind(xi - xi[1],yi - yi[1])
    vecf = cbind(xi - xi[length(xi)],yi - yi[length(yi)])
    
    leno = sqrt(sum(veco^2))
    leni = apply(veci,1,function(x) sqrt(x[1]^2 + x[2]^2))
    lenf = apply(vecf,1,function(x) sqrt(x[1]^2 + x[2]^2))
    
    theta = acos((leno^2 + leni^2 - lenf^2)/(2*leno*leni))
    
    eb = sin(theta)*leni
    eb[1] = 0
    
    ### Now, apply cushion to account for noisy elbow
    ebmx = max(eb,na.rm = TRUE)
    ebii = which(eb >= (ebmx*(1-cushion)))
    
    eb.idx = min(ebii)
    x.eb = xi[min(ebii)]
    
    if (verbose)
    {
      plot.new()
      par(mfrow = c(2,1))
      plot(xi,yi)
      plot(xi,eb);
      abline(v = x.eb,col = "red")
    }
    
    c("idx" = eb.idx,"x.elbow" = x.eb)
    
}


# get coordinates for centered labels 
get_module_label_coord <- function(coord,mods)
{
	do.call('rbind',lapply(mods,function(x,y) apply(y[match(x,rownames(y)),],2,function(q) median(q,na.rm = T)),y = coord))
}

embed_wt_network <- function(g,reduction = "umap",laplacian.type = "DAD",n.max = 50,get.elbow = TRUE)
{
	require(umap)
	require(Rtsne)
	cat("Calculate Laplacian and eigen-vectors....\n")
	laplacian = embed_laplacian_matrix(graph = g, no = n.max, weights = NULL, which = "lm",
	type = laplacian.type, scaled = TRUE)
	
	# find the elbow in the eigenvalues
	if (get.elbow)
	{
		cat("Find elbow in eigen-value plot and use dimension reduction...\n")
		nf = find_elbow_dist(x = 1:n.max,y = laplacian$D)
		nf = nf[1]
	}else{
		nf = min(c(ncol(laplacian$X),n.max))
	}
	cat(paste0("- Final dimension kept:",nf,"\n"))
	if (reduction == "umap")
	{
		out = umap::umap(laplacian$X[,1:nf])
		xy = out$layout;colnames(xy) = c("x","y")
		rownames(xy) = V(g)$name
	}
	
	if (reduction == "tsne")
	{
		out = Rtsne(X = laplacian$X[,1:nf],check_duplicates = FALSE)
		xy = out$Y;
		colnames(xy) = c("x","y")
		rownames(xy) = V(g)$name
	}
	
	output = list(coord = xy,laplacian.method = laplacian.type,nf = nf,eigen.vectors = laplacian$X,eigen.values = laplacian$D)
	return(output)
}


