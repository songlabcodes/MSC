evaluate_true_cluster_connectivity <- function(g,cls.t)
{
    require(igraph)
    require(Matrix)
    adj= get.adjacency(g)
    
    mem.t = membership.to.mem(cls.t)
    mem.t = mem.t[match(rownames(adj),rownames(mem.t)),]
    
    # cluster to cluster link ratio
    adj.count = t(mem.t) %*% adj %*% mem.t
    diag(adj.count) = 1/2 * diag(adj.count)
    
    # per group ratio
    inter.link = colSums(adj.count) - diag(adj.count)
    intra.link = diag(adj.count)
    
    # per node ratio
    adj.per.node = t(mem.t) %*% adj 
    intra.per.node = colSums(t(mem.t) * adj.per.node)
    inter.per.node = colSums(adj.per.node) - intra.per.node
    ratio.per.node = as.vector(t(mem.t) %*% cbind(intra.per.node/(inter.per.node + intra.per.node))/colSums(mem.t))
    names(ratio.per.node) = colnames(mem.t)
    # put together data
    grp.size = colSums(mem.t)
    link.df = data.frame(group.id = names(grp.size),group.size = grp.size,
                         intra.link = intra.link[names(grp.size)],inter.link = inter.link[names(grp.size)],
                         intra.inter.ratio = intra.link[names(grp.size)]/(inter.link[names(grp.size)] + intra.link[names(grp.size)]),
                         intra.inter.ratio.per.node = ratio.per.node[names(grp.size)]
    )
    return(link.df)
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

inclusion_rate <- function(mod.o,mod.f)
{
  # mod.o = gold standard
  # mod.f = inferred clusters
  lst = do.call('rbind',lapply(mod.f,function(x) sapply(mod.o,function(y,z) length(intersect(y,z))/length(z),z=x)))
  ir.vec = apply(lst,1,function(x) max(x,na.rm = T))
  mlen = sapply(mod.f,length)
  sum(ir.vec*mlen)/sum(mlen)
}

coverage_rate <- function(mod.o,mod.f)
{
  # mod.o = gold standard
  # mod.f = inferred clusters
  lst = do.call('rbind',lapply(mod.o,function(x) sapply(mod.f,function(y,z) length(intersect(y,z))/length(z),z=x)))
  cr.vec = apply(lst,1,function(x) max(x,na.rm = T))
  mlen = sapply(mod.o,length)
  sum(cr.vec*mlen)/sum(mlen)
}

accuracy_rate <- function(mod.o,mod.f)
{
  jac = calculate_jaccard(mod.o,mod.f)
  ac.vec = apply(jac,1,max)
  mlen = sapply(mod.o,length)
  sum(ac.vec*mlen)/sum(mlen)
}

accuracy_vector <- function(mod.o,mod.f)
{
  jac = calculate_jaccard(mod.o,mod.f)
  apply(jac,1,max)
  
}

evaluate_methods <- function(ct.lst,mod.lst)
{
  require(ggplot2)
  require(patchwork)
  require(cowplot)
  incl.rate <- covr.rate <- accr.rate <- rep(NA,length(mod.lst))
  names(incl.rate) <- names(covr.rate) <- names(accr.rate) <- names(mod.lst)
  for (i in 1:length(mod.lst))
  {
    incl.rate[i] = inclusion_rate(mod.o = ct.lst,mod.f = mod.lst[[i]])
    covr.rate[i] = coverage_rate(mod.o = ct.lst,mod.f = mod.lst[[i]])
    accr.rate[i] = accuracy_rate(mod.o = ct.lst,mod.f = mod.lst[[i]])
  }
  pdat = data.frame(method = names(mod.lst),incl.rate = incl.rate,covr.rate = covr.rate,accr.rate = accr.rate)
  p1 = ggplot(data = pdat,aes(x = method,y = incl.rate,fill = method)) + geom_bar(stat = "identity") + 
    labs(y = "Inclusion Rate") + 
    scale_x_discrete(limits = names(mod.lst)) + 
    guides(fill = "none") + 
    theme_classic() + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text(size = 17),
          axis.text.y = element_text(size = 14))
  p2 = ggplot(data = pdat,aes(x = method,y = covr.rate,fill = method)) + geom_bar(stat = "identity") + 
    labs(y = "Coverage Rate") + 
    scale_x_discrete(limits = names(mod.lst)) + 
    guides(fill = "none") + 
    scale_y_continuous(limits = c(0,1)) + 
    theme_classic() + 
    theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text(size = 17),
          axis.text.y = element_text(size = 14))
  p3 = ggplot(data = pdat,aes(x = method,y = accr.rate,fill = method)) + geom_bar(stat = "identity") + 
    labs(y = "Detection Accuracy") + 
    scale_x_discrete(limits = names(mod.lst)) + 
    guides(fill = "none") + 
    scale_y_continuous(limits = c(0,1)) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 14),axis.title.x = element_blank(),axis.title.y = element_text(size = 17),
          axis.text.y = element_text(size = 14))
  
  library(patchwork)
  perf.plot = p1/p2/p3
  return(list(pobj = perf.plot,pdata = pdat))
}
