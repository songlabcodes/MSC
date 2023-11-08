##### Clustering functions
# functions
membership.to.mem <- function(membership)
  {
    cls <- setdiff(unique(membership),NA)
    mem <- Matrix::Matrix(0,nrow = length(membership),ncol = length(cls))
    for (i in 1:length(cls)) mem[which(membership == cls[i]),i] <- 1;
    colnames(mem) <- cls;
    rownames(mem) <- names(membership)
    return(mem)
  }
  
  ### compactness functions
  # shortest path calculation
  get.sp <- function(g,d.func = function(x){1-E(x)$weight})
  {
    dg <- g;
    E(dg)$weight <- d.func(g);
    sp.dist <- distances(graph = dg, v = V(dg), to = V(dg), mode = "all", weights = NULL, algorithm = "bellman-ford")
    return(sp.dist)
  }
  
  # compactness score
  compact.indiv <- function(v,D,alpha)
  {
    if (length(v) > 1)
    {
      d <- D[match(v,rownames(D)),match(v,colnames(D))]
      SPD.mu <- median(d[upper.tri(d)],na.rm = T)
      out <- SPD.mu/(log(length(v))^alpha)
    }else{
      out <- Inf;
    }
    return(out)
  }
  
  compact.indiv.graph <- function(g,d.func,alpha)
  {
    out = Inf
    if (vcount(g) > 1)
    {
      SPD.mu = mean_distance(graph=g,weights = d.func(E(g)$weight))
      out = SPD.mu/(log(vcount(g))^alpha)
    }
    return(out)
  }
  
  find_optimal_solution.leiden.v2 <- function(g,resol.vec = seq(0.1,2,0.0125),pcut= 0.05,n.perm = 10,min.width = 0.2,seed = 1234)
  {
    
    grand.mean <- function(M, N) {weighted.mean(M, N)}
    grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                            weighted.mean(M, N)^2)}
    # m1, m2: the sample means
    # s1, s2: the sample standard deviations
    # n1, n2: the same sizes
    # m0: the null value for the difference in means to be tested for. Default is 0. 
    # equal.variance: whether or not to assume equal variance. Default is FALSE. 
    t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
    {
      if( equal.variance==FALSE ) 
      {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
      } else
      {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
        df <- n1+n2-2
      }      
      t <- (m1-m2-m0)/se 
      dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
      names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
      return(dat) 
    }
    
    output = NULL 
    clsm = kratio = matrix(0,nrow = vcount(g),ncol = 0);
    rownames(clsm) = rownames(kratio) = V(g)$name
    cls.f = NULL
    adj = as_adj(graph = g,attr = "weight")
    kr=vector("list",length(resol.vec))
    tstat.mat = matrix(nrow = 0,ncol = 5);
    qval = rep(NA,length(resol.vec))
    set.seed(seed)
    for (ri in 1:length(resol.vec))
    {
      cout = cluster_leiden(graph = g,objective_function = "modularity",resolution_parameter = resol.vec[ri])
      qval[ri] = cout$quality
      clsm = cbind(clsm,cout$membership)
      
      mem = membership.to.mem(membership = cout$membership)
      k.global <- Matrix::rowSums(adj,na.rm = T)
      mm = (Matrix::t(mem) %*% adj) * Matrix::t(mem)
      #print(dim(mm))
      k.in <- Matrix::colSums(mm)
      so = (k.global - k.in)/k.global;
      
    
      sr = lapply(1:n.perm,function(n) {
        mem = membership.to.mem(membership = sample(cout$membership,size = length(cout$membership)))
        k.global <- Matrix::rowSums(adj,na.rm = T)
        k.in <- Matrix::colSums((Matrix::t(mem) %*% adj) * Matrix::t(mem))
        (k.global - k.in)/k.global;
      })
      
      df.r = data.frame(n = sapply(sr,length),mean = sapply(sr,mean),sd = sapply(sr,sd))
      
      mu.r = grand.mean(M = df.r$mean,N = df.r$n)
      sd.r = grand.sd(S = df.r$sd,M = df.r$mean,N = df.r$n)
      
      mu.o = mean(so);sd.o = sd(so)
      
      t.stat = t.test2(m1 = mu.o,m2 = mu.r,s1 = sd.o,s2 = sd.r,n1 = nrow(mem),n2 = nrow(mem),m0=0,equal.variance=FALSE)
      t.stat = c(t.stat,"fold.change" = mu.o/mu.r)
      tstat.mat = rbind(tstat.mat,t.stat)
      rm(t.stat)
    }
    
    # check
    if (any(is.nan(tstat.mat[,3]) | qval < 0)) 
    {
      ii = !is.nan(tstat.mat[,3]) & qval > 0 & tstat.mat[,4] < pcut
      tstat.mat = rbind(tstat.mat[ii,])
      resol.vec = resol.vec[ii]
    }
    
    output = NULL
    if (nrow(tstat.mat) > 5)
    {
      do.check = FALSE
      if (length(unique(tstat.mat[,5])) > 1)
      {
        # check if the t-statistics do have steps in the curves: if there is a step, the distribution shouldn't follow normal distribution
        step.check = shapiro.test(x = tstat.mat[,5])
        if (step.check$p.value < pcut) do.check = TRUE
        
      }else{
        ii = ceiling(length(resol.vec)/2)
        resol.f = resol.vec[ii]
        cls.f = clsm[,ii]
        output = list(optimal.modules = factor(cls.f),cluster.matrix = clsm,
                      resolution.tested = resol.vec,
                      resolution.used = resol.f)
      }
      if (do.check)
      {
        cat("Step detected...\n")
        tree <- rpart(-tstat.mat[,5] ~ resol.vec)
        
        # find the break points leading to biggest gap
        yfit = predict(tree, data.frame(x=resol.vec))
        
        if (any(abs(diff(yfit)) > 1E-320))
        {
          break.points = which(abs(diff(yfit)) > 1E-320)
          break.dir = sign(diff(yfit))[break.points]
          
          init = 1
          break.points = c(break.points,length(yfit))
          intervals = matrix(0,nrow = 0,ncol = 2)
          for (i in 1:length(break.points))
          {
            intervals = rbind(intervals,c(init,break.points[i]))
            init = break.points[i] + 1
          }
          rm(init)
          interval.val = yfit[intervals[,1]]
          
          interval.data = data.frame(xo = intervals[,1],xf = intervals[,2],height = interval.val,
                                     break.dir = c(0,break.dir))
          nc = rep(NA,nrow(interval.data))
          for (i in 1:nrow(interval.data)) nc[i] = mean(apply(clsm[,interval.data$xo[i]:interval.data$xf[i]],2,function(x) length(unique(x))))
          interval.data$nc = nc
          
          # clear into monotonic descending data
          init = interval.data$height[1]
          bad.row = c()
          for (i in 2:nrow(interval.data))
          {
            if (interval.data$height[i] > init) bad.row = c(bad.row,i)
            if (interval.data$height[i] < init) init = interval.data$height[i]
          }
          vec = rep(TRUE,nrow(interval.data));vec[bad.row] = FALSE
          interval.data$is.descending = vec
          
          vec = rep(FALSE,nrow(interval.data))
          interval.data$resol.width = resol.vec[interval.data$xf] - resol.vec[interval.data$xo]
          
          # get cluster number call
          interval.data$n.cluster = apply(do.call('cbind',interval.data[,1:2]),1,function(x) median(apply(clsm[,x[1]:x[2]],2,function(y) length(unique(y)))))
          
          # get the initial drop
          vi = min(which(interval.data$is.descending & interval.data$resol.width >= min.width & interval.data$n.cluster > 1))
          if (!is.infinite(vi)) 
          {
            resol.f = resol.vec[median(interval.data$xo[vi]:interval.data$xf[vi])]
            cls.f = clsm[,which(resol.vec == resol.f)]
          }else{
            resol.f = NULL
          }
          
        }else{
          resol.f = NULL
        }
        
        if (!is.null(resol.f))
        {
          # plot out fitted steps
          plot.new()
          plot(resol.vec,-tstat.mat[,5],xlab = "Resolution (\u03B3)",ylab = "")
          lines(resol.vec,predict(tree, data.frame(x=resol.vec)),col = "blue")
          abline(v = resol.f,col = "red",lty = 2)
          
          cls.f = clsm[,which(resol.vec == resol.f)]
          output = list(optimal.modules = factor(cls.f),cluster.matrix = clsm,
                        resolution.tested = resol.vec,
                        resolution.used = resol.f)
          
        }
        
        
      }else{
        ii = ceiling(length(resol.vec)/2)
        resol.f = resol.vec[ii]
        cls.f = clsm[,ii]
        output = list(optimal.modules = factor(cls.f),cluster.matrix = clsm,
                      resolution.tested = resol.vec,
                      resolution.used = resol.f)
      }
    }
    return(output)
  }
  
  #' @title Perform adaptive split of a network (AdaptSplit)
  #' @name split_network
  #' @docType package
  #' @description Sweeps over the clustering resolution in \code{[0.1,2]} to find the most granular clustering results. 
  #' @param g An igraph object containing the network
  #' @param mi An integer. The clusters will be indexed, starting from the integer specified by mi.  
  #' @param d.func a function to convert the edge weights into edge function. 
  #' @param alpha Compactness resolution parameter. This value is used to calculate the compactness of the child clusters obtained from the split. 
  #' @param node.perturb.prop Proportion of nodes to randomly shuffle to calculate the reference distribution of intra-cluster connectivity. 
  #' @return Returns a list containing clustering results with optimal.modules (the final results from adapt split), resolution.tested (clustering resolution sweeped), resolution.used (the final value of clustering resolution), and module.table containing the compactness and intra-cluster connectivity stats.
  #' @examples 
  #' data(pbmc_8k_msc_results)
  #' split_network(g = pbmc_8k_msc_results$cell.network,mi = 0)
  #' @export
  NULL
  split_network <- function(g,mi,d.func = function(x) x,alpha = 1,node.perturb.prop = 0.05)
  {
    require(rpart)
    find_optimal_solution.leiden.v2 <- function(g,resol.vec = seq(0.1,2,0.0125),pcut= 0.05,n.perm = 10,min.width = 0.2,seed = 1234)
    {
      
      grand.mean <- function(M, N) {weighted.mean(M, N)}
      grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                              weighted.mean(M, N)^2)}
      # m1, m2: the sample means
      # s1, s2: the sample standard deviations
      # n1, n2: the same sizes
      # m0: the null value for the difference in means to be tested for. Default is 0. 
      # equal.variance: whether or not to assume equal variance. Default is FALSE. 
      t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
      {
        if( equal.variance==FALSE ) 
        {
          se <- sqrt( (s1^2/n1) + (s2^2/n2) )
          # welch-satterthwaite df
          df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
        } else
        {
          # pooled standard deviation, scaled by the sample sizes
          se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
          df <- n1+n2-2
        }      
        t <- (m1-m2-m0)/se 
        dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
        names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
        return(dat) 
      }
      
      output = NULL 
      clsm = kratio = matrix(0,nrow = vcount(g),ncol = 0);
      rownames(clsm) = rownames(kratio) = V(g)$name
      cls.f = NULL
      adj = as_adj(graph = g,attr = "weight")
      kr=vector("list",length(resol.vec))
      tstat.mat = matrix(nrow = 0,ncol = 5);
      qval = rep(NA,length(resol.vec))
      set.seed(seed)
      for (ri in 1:length(resol.vec))
      {
        cout = cluster_leiden(graph = g,objective_function = "modularity",resolution_parameter = resol.vec[ri])
        qval[ri] = cout$quality
        clsm = cbind(clsm,cout$membership)
        
        mem = membership.to.mem(membership = cout$membership)
        k.global <- Matrix::rowSums(adj,na.rm = T)
        mm = (Matrix::t(mem) %*% adj) * Matrix::t(mem)
        #print(dim(mm))
        k.in <- Matrix::colSums(mm)
        so = (k.global - k.in)/k.global;
        
        
        sr = lapply(1:n.perm,function(n) {
          mem = membership.to.mem(membership = sample(cout$membership,size = length(cout$membership)))
          k.global <- Matrix::rowSums(adj,na.rm = T)
          k.in <- Matrix::colSums((Matrix::t(mem) %*% adj) * Matrix::t(mem))
          (k.global - k.in)/k.global;
        })
        
        df.r = data.frame(n = sapply(sr,length),mean = sapply(sr,mean),sd = sapply(sr,sd))
        
        mu.r = grand.mean(M = df.r$mean,N = df.r$n)
        sd.r = grand.sd(S = df.r$sd,M = df.r$mean,N = df.r$n)
        
        mu.o = mean(so);sd.o = sd(so)
        
        t.stat = t.test2(m1 = mu.o,m2 = mu.r,s1 = sd.o,s2 = sd.r,n1 = nrow(mem),n2 = nrow(mem),m0=0,equal.variance=FALSE)
        t.stat = c(t.stat,"fold.change" = mu.o/mu.r)
        tstat.mat = rbind(tstat.mat,t.stat)
        rm(t.stat)
      }
      
      # check
      if (any(is.nan(tstat.mat[,3]) | qval < 0)) 
      {
        ii = !is.nan(tstat.mat[,3]) & qval > 0 & tstat.mat[,4] < pcut
        tstat.mat = rbind(tstat.mat[ii,])
        resol.vec = resol.vec[ii]
      }
      
      output = NULL
      if (nrow(tstat.mat) > 5)
      {
        do.check = FALSE
        if (length(unique(tstat.mat[,5])) > 1)
        {
          # check if the t-statistics do have steps in the curves: if there is a step, the distribution shouldn't follow normal distribution
          step.check = shapiro.test(x = tstat.mat[,5])
          if (step.check$p.value < pcut) do.check = TRUE
          
        }else{
          ii = ceiling(length(resol.vec)/2)
          resol.f = resol.vec[ii]
          cls.f = clsm[,ii]
          output = list(optimal.modules = factor(cls.f),cluster.matrix = clsm,
                        resolution.tested = resol.vec,
                        resolution.used = resol.f)
        }
        if (do.check)
        {
          cat("Step detected...\n")
          tree <- rpart(-tstat.mat[,5] ~ resol.vec)
          
          # find the break points leading to biggest gap
          yfit = predict(tree, data.frame(x=resol.vec))
          
          if (any(abs(diff(yfit)) > 1E-320))
          {
            break.points = which(abs(diff(yfit)) > 1E-320)
            break.dir = sign(diff(yfit))[break.points]
            
            init = 1
            break.points = c(break.points,length(yfit))
            intervals = matrix(0,nrow = 0,ncol = 2)
            for (i in 1:length(break.points))
            {
              intervals = rbind(intervals,c(init,break.points[i]))
              init = break.points[i] + 1
            }
            rm(init)
            interval.val = yfit[intervals[,1]]
            
            interval.data = data.frame(xo = intervals[,1],xf = intervals[,2],height = interval.val,
                                       break.dir = c(0,break.dir))
            nc = rep(NA,nrow(interval.data))
            for (i in 1:nrow(interval.data)) nc[i] = mean(apply(clsm[,interval.data$xo[i]:interval.data$xf[i]],2,function(x) length(unique(x))))
            interval.data$nc = nc
            
            # clear into monotonic descending data
            init = interval.data$height[1]
            bad.row = c()
            for (i in 2:nrow(interval.data))
            {
              if (interval.data$height[i] > init) bad.row = c(bad.row,i)
              if (interval.data$height[i] < init) init = interval.data$height[i]
            }
            vec = rep(TRUE,nrow(interval.data));vec[bad.row] = FALSE
            interval.data$is.descending = vec
            
            vec = rep(FALSE,nrow(interval.data))
            interval.data$resol.width = resol.vec[interval.data$xf] - resol.vec[interval.data$xo]
            
            # get cluster number call
            interval.data$n.cluster = apply(do.call('cbind',interval.data[,1:2]),1,function(x) median(apply(clsm[,x[1]:x[2]],2,function(y) length(unique(y)))))
            
            # get the initial drop
            vi = min(which(interval.data$is.descending & interval.data$resol.width >= min.width & interval.data$n.cluster > 1))
            if (!is.infinite(vi)) 
            {
              resol.f = resol.vec[median(interval.data$xo[vi]:interval.data$xf[vi])]
              cls.f = clsm[,which(resol.vec == resol.f)]
            }else{
              resol.f = NULL
            }
            
          }else{
            resol.f = NULL
          }
          
          if (!is.null(resol.f))
          {
            # plot out fitted steps
            plot.new()
            plot(resol.vec,-tstat.mat[,5],xlab = "Resolution (\u03B3)",ylab = "")
            lines(resol.vec,predict(tree, data.frame(x=resol.vec)),col = "blue")
            abline(v = resol.f,col = "red",lty = 2)
            
            cls.f = clsm[,which(resol.vec == resol.f)]
            output = list(optimal.modules = factor(cls.f),cluster.matrix = clsm,
                          resolution.tested = resol.vec,
                          resolution.used = resol.f)
            
          }
          
          
        }else{
          ii = ceiling(length(resol.vec)/2)
          resol.f = resol.vec[ii]
          cls.f = clsm[,ii]
          output = list(optimal.modules = factor(cls.f),cluster.matrix = clsm,
                        resolution.tested = resol.vec,
                        resolution.used = resol.f)
        }
      }
      return(output)
    }
    membership.to.mem <- function(membership)
    {
      cls <- setdiff(unique(membership),NA)
      mem <- Matrix::Matrix(0,nrow = length(membership),ncol = length(cls))
      for (i in 1:length(cls)) mem[which(membership == cls[i]),i] <- 1;
      colnames(mem) <- cls;
      rownames(mem) <- names(membership)
      return(mem)
    }
    convert_mat_to_edgelist <- function(m)
    {
      ## quickly convert matrix to three column list
      ij = do.call('rbind',lapply(1:(ncol(m)-1),function(x,n) {rr = rep(x);cc = (x+1):n;cbind(rr,cc)},n = ncol(m)));
      colnames(ij) = c("ri","ci")
      w = m[ij_to_lin(rr = ij[,1],cc = ij[,2],n = ncol(m))]
      cbind(ij,w)
    }
    cat("# Commence coarse-grained structure detection...\n")
    output = NULL
    #g = subgraph(graph = g,vids = mods$M12)
    
    set.seed(1234)
    
    clus.res = find_optimal_solution.leiden.v2(g = g,min.width = 0.1)
    
    if (!is.null(clus.res))
    {
      clus.resol = clus.res$resolution.used
      clsvec = factor(clus.res$optimal.modules);
      
      #### Now, check out what criteria work for terminating cluster search
      
      # update module names with M prefix
      levels(clsvec) = paste0("M",(mi+1):(mi+length(levels(clsvec))));
      clus.res$optimal.modules = clsvec
      mods = split(names(clsvec),clsvec)
      cat("- clusters:\n");
      print(sapply(mods,length))
      
      ###### compactness
      cat(paste0("- evaluate compactness with alpha = ",alpha,"\n"))
      
      # global mst for compactness architecture
      g.mst = mst(graph = g,weights = d.func(E(g)$weight))
      sp.dist = distances(graph = g.mst, v = V(g.mst), to = V(g.mst), mode = "all", weights = d.func(E(g.mst)$weight), algorithm = "bellman-ford")
      comp.o = mean(sp.dist[upper.tri(sp.dist)])/log(vcount(g.mst))^alpha
      
      # get module-wise compactness
      compact.m = sapply(mods,function(x,sp.dist,alpha) 
      {
        spd = sp.dist[match(x,rownames(sp.dist)),match(x,colnames(sp.dist))]
        mean(spd[upper.tri(spd)])/log(length(x))^alpha
      },sp.dist = sp.dist,alpha = alpha)
      
      # now, get random network 
      #compact.m = c(compact.m,comp.o)
      cat("-- Check random compactness for controls...\n")
      ns = 100
      pval=rand.compact = rep(NA,length(mods))
      for (i in 1:length(mods))
      {
        rand.comp= rep(NA,ns)
        for (n in 1:ns)
        {
          ii = sample(1:nrow(sp.dist),length(mods[[i]]))
          spd = sp.dist[ii,ii];
          rand.comp[n] = mean(spd[upper.tri(spd)])/log(length(mods[[i]]))^alpha
          
        }
        pval[i] =sum(rand.comp < compact.m[i])/ns
        rand.compact[i] = mean(rand.comp)
      }
      
      ### edge density
      cat("- evaluate intra-cluster connectivity\n")
      calculate_intra_link_ratio <- function(g,clsvec)
      {
        adjm = get.adjacency(g)
        mem = membership.to.mem(clsvec)
        adjm = adjm[match(rownames(mem),rownames(adjm)),match(rownames(mem),colnames(adjm))]
        link.mat = Matrix::t(mem) %*% adjm %*% mem;
        link.mat = as.matrix(link.mat)
        
        diag(link.mat) = diag(link.mat)/2
        intra.link.ratio = diag(link.mat)/rowSums(link.mat)
        return(intra.link.ratio)
      }
      
      intra.link.ratio = calculate_intra_link_ratio(g,clsvec = clus.res$optimal.modules)
      #intra.link = sapply(mods,function(x,g) ecount(induced.subgraph(graph = g,vids = x)),g = g)
      #intra.link.ratio = intra.link/ecount(g)
      
      cat("-- Check random connectivity for controls...\n")
      ns = 100
      intra.link.ratio.rand = matrix(0,nrow = ns,ncol = length(mods))
      pval.link = rep(NA,length(mods))
      for (n in 1:ns)
      {
        clus.rand= clus.res$optimal.modules;
        #names(clus.rand) = sample(names(clus.rand),length(clus.rand))
        # shuffle only 10% for subtle changes in the random
        ii = sample(1:length(clus.rand),ceiling(length(clus.rand)*node.perturb.prop))
        names(clus.rand)[ii] = sample(names(clus.rand)[ii],length(ii))
        #intra.link.rand = sapply(split(names(clus.rand),clus.rand),function(x,g) ecount(induced.subgraph(graph = g,vids = x)),g = g)
        intra.link.ratio.rand[n,] = calculate_intra_link_ratio(g,clsvec = clus.rand)
      }
      for (i in 1:length(pval.link)) pval.link[i] = sum(intra.link.ratio.rand[,i] > intra.link.ratio[i])

      ### get data
      cat("# Collect results\n")
      module.stat = data.frame(cluster.name = names(mods),module.size = sapply(mods,length),cluster.resolution = rep(clus.res$resolution.used,length(mods)),
                               compactness.resolution = rep(alpha,length(mods)),
                               compactness.parent = rep(comp.o,length(mods)),
                               compactness.module = compact.m,compactness.module.random= rand.compact,compactness.pvalue = pval,is.compact = compact.m < comp.o & pval < 0.05,
                               intra.link.ratio = intra.link.ratio,intra.link.ratio.random = apply(intra.link.ratio.rand,2,mean),intra.link.pvalue = pval.link,is.linked = pval.link < 0.05
      )
      
      output = c(clus.res,list(module.table = module.stat))
    }
    
    return(output)
  }
  
  #' @title Calculate the alpha value for compactness 
  #' @name identify_scale_parameter
  #' @docType package
  #' @description Estimates the compactness parameter alpha via subnetwork sampling. 
  #' @param g An igraph object containing the network
  #' @param nv An integer. Number of random nodes to sample for neighborhood search.
  #' @param nr.max Neighborhood layer to expand.
  #' @param d.func a function to convert the edge weights into edge function. 
  #' @param mode A character value from either of "diameter" or "mean". Decide to estimate the compactness based on network diameter (i.e. the longest of the pairwise shortest path distance), or mean (i.e. mean of all pairwise shortest path distance). 
  #' @return Returns a numeric value for the alpha value. 
  #' @examples 
  #' data(pbmc_8k_msc_results)
  #' identify_scale_parameter(pbmc_8k_msc_results$cell.network,mode = "diameter")
  #' @export
  NULL
  identify_scale_parameter <- function(g,nv = 100,nr.max = 3,seed = 1234,d.func = function(x) x,mode = c("diameter","mean"))
  {
    overall.size = alpha.val = c()
    set.seed(seed)
    for (nr in 1:nr.max)
    {
      vo = sample(V(g)$name,min(c(nv,vcount(g))))
      
      nlst = make_ego_graph(g,order = nr,nodes = vo,mode = "all")
      nlst = lapply(nlst,function(x,d.func) mst(graph = x,weights = d.func(E(x)$weight)),d.func = d.func)
      if (mode == "mean") diam = sapply(nlst,function(x) mean_distance(graph = x,directed = FALSE))
      if (mode == "diameter") diam = sapply(nlst,function(x) diameter(graph = x,directed = FALSE))
      gsize = sapply(nlst,vcount)
      
      overall.size = c(overall.size, gsize);
      alpha.val = c(alpha.val,log(diam)/log(log(gsize)))
    }
    
    plot.new()
    plot(overall.size,alpha.val,log = "x",ylab = "Scale parameter (\u03B1)",xlab = "Module Size");
    
    if (sum(overall.size > 40) > 10)
    {
      ocut = 40
      abline(h = median(alpha.val[overall.size > 40]),col = "red")
      alpha.o = median(alpha.val[overall.size > 40])
    }else{
      ocut = quantile(overall.size,0.8)
      abline(h = median(alpha.val[overall.size >= ocut]),col = "red")
      alpha.o = median(alpha.val[overall.size >= ocut])
    }
    #alpha.o = quantile(alpha.val[overall.size > 20],0.05)
    cat("Scale parameter:\n")
    print(alpha.o)
    
    return(alpha.o)
    
  }
  
  #' @title Iterative top-down clustering workflow 
  #' @name iterative_clustering
  #' @docType package
  #' @description Given input cell network, performs the iterative clustering to probe the multi-scale structure. 
  #' @param g.in An igraph object containing the network
  #' @param min.size An integer. The minumum cluster size threshold.
  #' @param alpha The alpha exponent for compactness calculation. Can be obtained from identify_scale_parameter().
  #' @param d.func a function to convert the edge weights into edge function. Default is d(x) = sqrt(2*(1-x)), the metric distance for correlation coefficient. 
  #' @param pcut A numeric. The p-value to identify clusters with significantly high intra-cluster connectivity. Default is 0.05.  
  #' @param seed Seed value to initialize the random sampling. Default is 1234.
  #' @return Returns a numeric value for the alpha value. 
  #' @seealso [identify_scale_parameter()]
  #' @examples 
  #' \donttest{
  #' data(pbmc_8k_msc_results)
  #' alpha.val = identify_scale_parameter(pbmc_8k_msc_results$cell.network,mode = "diameter")
  #' iter.res = iterative_clustering(g.in = pbmc_8k_msc_results$cell.network,alpha = alpha.val)
  #' }
  #' @export
  iterative_clustering <- function(g.in,min.size = 10,d.func = function(x) sqrt(2*(1-x)),alpha = 1,pcut = 0.05,seed = 1234)
  {
    # initialize
    go = list(M0 = g.in)
    do.run = TRUE
    mi = 0
    n.iter = 0
    
    modules = list(M0 = V(g.in)$name)
    
    module.table = data.frame()
    set.seed(seed)
    while (do.run)
    {
      n.iter = n.iter + 1
      cat(paste0("##### Multiscale clustering iteration.",n.iter,"\n"))
      cat(paste0(" - current number of identified modules:",nrow(module.table),"\n"))
      cat(paste0(" - number of parent modules to split:",length(go),"\n"))
      cat(paste0(" - parent module names:",paste(names(go),collapse = ","),"\n"))
      
      gn = list()
      indiv.run = rep(TRUE,length(go))
      for (i in 1:length(go))
      {
        is.single.graph = FALSE
        gi = go[[i]]
        if (is.connected(gi))
        {
          is.single.graph = TRUE
          cout = components(gi)
          if (sum(cout$csize > log(vcount(gi))) == 1)
          {
            ii = which(cout$csize > log(vcount(gi)))
            gi = induced_subgraph(graph = gi,vids = names(cout$membership)[cout$membership == ii])
            
          }
          
        }
        
        # if connected graph, run split_network routine. If not, get connected components as the first layer of results
        out = NULL
        if (is.single.graph)
        {
          out = split_network(g = gi,mi = mi,d.func = d.func,alpha = alpha)
          
        }else{
          cout = components(gi);
          clsvec = factor(cout$membership)
          levels(clsvec) = paste0("M",(mi+1):(mi+length(levels(clsvec))));
          mods = split(names(clsvec),clsvec)
          
          # global mst for compactness architecture
          compact.m = sapply(mods,function(m) {
            g = induced_subgraph(gi,vids = m)
            g.mst = mst(graph = g,weights = d.func(E(g)$weight))
            sp.dist = distances(graph = g.mst, v = V(g.mst), to = V(g.mst), mode = "all", weights = d.func(E(g.mst)$weight), algorithm = "bellman-ford")
            mean(sp.dist[upper.tri(sp.dist)])/log(vcount(g.mst))^alpha
            
          })
          
          module.stat = data.frame(cluster.name = names(mods),module.size = sapply(mods,length),cluster.resolution = rep(NA,length(mods)),
                                   compactness.resolution = rep(alpha,length(mods)),
                                   compactness.parent = rep(Inf,length(mods)),
                                   compactness.module = compact.m,compactness.module.random= rep(Inf,length(mods)),compactness.pvalue = rep(0,length(mods)),is.compact = rep(TRUE,length(mods)),
                                   intra.link.ratio = rep(1,length(mods)),intra.link.ratio.random = rep(NA,length(mods)),intra.link.pvalue = rep(0,length(mods)),is.linked = rep(TRUE,length(mods)))
          
          out = list(optimal.modules = clsvec,module.table = module.stat)
        }
        
        if (!is.null(out) & length(unique(out$optimal.modules)) > 1)
        {
          mods = split(names(out$optimal.modules),factor(out$optimal.modules))
          mi = mi + nrow(out$module.table)
          
          # check if smaller than infimum compactness across all parents
          tbl = out$module.table;
          tbl$parent.cluster = rep(names(go)[i],nrow(tbl));
          cols = c("cluster.name","module.size","cluster.resolution","compactness.resolution","parent.cluster","compactness.parent","compactness.module",
                   "compactness.pvalue","intra.link.ratio","intra.link.pvalue")
          
          if (nrow(module.table) > 0)
          {
            mtbl = rbind.data.frame(module.table[,cols],tbl[,cols])
          }else{
            mtbl = rbind.data.frame(tbl[,cols])
          }
          
          comp.o = subset(mtbl,parent.cluster == "M0")$compactness.parent[1]
          gh = graph.data.frame(mtbl[,c("cluster.name","parent.cluster")],directed = TRUE)
          prnts = lapply(mtbl$cluster.name,function(x) setdiff(names(ego(graph = gh,order = 5,nodes = x,mode = "out")[[1]]),x))
          prnts.comp = sapply(prnts,function(x) {
            out = c(comp.o,subset(mtbl,cluster.name %in% x)$compactness.module)
            val = min(out,na.rm = T)
            return(val)
          })
          if (any(is.infinite(prnts.comp))) prnts.comp[is.infinite(prnts.comp)] = subset(mtbl,parent.cluster == "M0")$compactness.parent
          mtbl$compactness.parent.infimum = prnts.comp
          mtbl$is.compact = mtbl$compactness.module < mtbl$compactness.parent.infimum & mtbl$compactness.pvalue < pcut
          mtbl$is.coherent = mtbl$intra.link.pvalue < pcut
          mtbl$is.valid = (mtbl$is.compact | mtbl$is.coherent) & mtbl$module.size >= min.size
          mtbl = subset(mtbl,cluster.name %in% names(mods))
          if (any(mtbl$is.valid))
          {
            # update modules
            modules = c(modules,mods)
            module.table = rbind.data.frame(module.table,mtbl)
            
            gc = lapply(mods[subset(mtbl,is.valid)$cluster.name],function(x) subgraph(go[[i]],vids = x))
            gn = c(gn,gc)
          }else{
            indiv.run[i] = FALSE
          }
        }else{
          indiv.run[i] = FALSE
        }
        
        
      }
      
      gn = gn[sapply(gn,vcount) >= min.size]
      if (any(indiv.run) & length(gn) > 0)
      {
        do.run = TRUE
        go = gn
      }else{
        do.run = FALSE
      }
      
    }
    
    output = list(modules = modules,module.table = module.table,cell.network = g.in,alpha.value = alpha)
    return(output)
  }
  
  #' @title Parallelized, iterative top-down clustering workflow 
  #' @name iterative_clustering.par
  #' @docType package
  #' @description Given input cell network, performs the iterative clustering to probe the multi-scale structure with parallelization. 
  #' @param g.in An igraph object containing the network
  #' @param min.size An integer. The minumum cluster size threshold.Default is 10. 
  #' @param alpha The alpha exponent for compactness calculation. Can be obtained from identify_scale_parameter().
  #' @param d.func a function to convert the edge weights into edge function. Default is d(x) = sqrt(2*(1-x)), the metric distance for correlation coefficient. 
  #' @param pcut A numeric. The p-value to identify clusters with significantly high intra-cluster connectivity. Default is 0.05.  
  #' @param seed Seed value to initialize the random sampling. Default is 1234.
  #' @param valid.gate Filtering method to identify significant clusters. "density" compares intra-cluster link density between parent and child clusters. "compact" compared compactness with the parent clusters. "both" uses "density" and "compact" filters simultaneously. Default is "both".
  #' @param n.cores An integer for the number of cores for parallelization 
  #' @return Returns a numeric value for the alpha value. 
  #' @seealso [identify_scale_parameter()]
  #' @examples 
  #' \donttest{
  #' data(pbmc_8k_msc_results)
  #' alpha.val = identify_scale_parameter(pbmc_8k_msc_results$cell.network,mode = "diameter")
  #' iter.res = iterative_clustering.par(g.in = pbmc_8k_msc_results$cell.network,alpha = alpha.val,n.cores = 4,valid.gate = "both")
  #' }
  #' @export
  iterative_clustering.par <- function(g.in,min.size = 10,
                                   d.func = function(x) sqrt(2*(1-x)),
                                   alpha = 1,pcut = 0.05,seed = 1234,
                                   valid.gate = c("density","compact","both"),
                                   n.cores = 4)
  {
    require(doParallel)
    # register parallel backend
    if (getDoParWorkers() == 1 & n.cores > 1)
    {
      cl <- parallel::makeCluster(n.cores)
      registerDoParallel(cl)
      # check how many workers are there
      cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
    }
    
    # initialize
    go = list(M0 = g.in)
    do.run = TRUE
    mi = 0
    n.iter = 0
    
    modules = list(M0 = V(g.in)$name)
    
    module.table = data.frame()
    set.seed(seed)
    while (do.run)
    {
      n.iter = n.iter + 1
      cat(paste0("##### Multiscale clustering iteration.",n.iter,"\n"))
      cat(paste0(" - current number of identified modules:",nrow(module.table),"\n"))
      cat(paste0(" - number of parent modules to split:",length(go),"\n"))
      cat(paste0(" - parent module names:",paste(names(go),collapse = ","),"\n"))
      
      gn = list()
      indiv.run = rep(TRUE,length(go))
      
      # run clustering in parallel
      out.lst <- foreach(i=1:length(go),.packages = c("igraph","MEGENA","Matrix","rpart"),.export = c("convert_mat_to_edgelist","find_optimal_solution.leiden.v2","split_network","membership.to.mem")) %dopar% {
        require(rpart)
        is.single.graph = FALSE
        gi = go[[i]]
        if (is.connected(gi))
        {
          is.single.graph = TRUE
          cout = components(gi)
          if (sum(cout$csize > log(vcount(gi))) == 1)
          {
            ii = which(cout$csize > log(vcount(gi)))
            gi = induced_subgraph(graph = gi,vids = names(cout$membership)[cout$membership == ii])
            
          }
          
        }
        
        # if connected graph, run split_network routine. If not, get connected components as the first layer of results
        out = NULL
        if (is.single.graph)
        {
          out = split_network(g = gi,mi = mi,d.func = d.func,alpha = alpha)
          
        }else{
          cout = components(gi);
          clsvec = factor(cout$membership)
          levels(clsvec) = paste0("M",(mi+1):(mi+length(levels(clsvec))));
          mods = split(names(clsvec),clsvec)
          
          # global mst for compactness architecture
          compact.m = sapply(mods,function(m) {
            g = induced_subgraph(gi,vids = m)
            g.mst = mst(graph = g,weights = d.func(E(g)$weight))
            sp.dist = distances(graph = g.mst, v = V(g.mst), to = V(g.mst), mode = "all", weights = d.func(E(g.mst)$weight), algorithm = "bellman-ford")
            mean(sp.dist[upper.tri(sp.dist)])/log(vcount(g.mst))^alpha
            
          })
          
          module.stat = data.frame(cluster.name = names(mods),module.size = sapply(mods,length),cluster.resolution = rep(NA,length(mods)),
                                   compactness.resolution = rep(alpha,length(mods)),
                                   compactness.parent = rep(Inf,length(mods)),
                                   compactness.module = compact.m,compactness.module.random= rep(Inf,length(mods)),compactness.pvalue = rep(0,length(mods)),is.compact = rep(TRUE,length(mods)),
                                   intra.link.ratio = rep(1,length(mods)),intra.link.ratio.random = rep(NA,length(mods)),intra.link.pvalue = rep(0,length(mods)),is.linked = rep(TRUE,length(mods)))
                                   
          out = list(optimal.modules = clsvec,module.table = module.stat)
        }
        return(out)
      }
      
      for (i in 1:length(out.lst))
      {
        out = out.lst[[i]]
        if (!is.null(out) & length(unique(out$optimal.modules)) > 1)
        {
          mods = split(names(out$optimal.modules),factor(out$optimal.modules))
          #mi = mi + nrow(out$module.table)
          
          # check if smaller than infimum compactness across all parents
          tbl = out$module.table;
          tbl$parent.cluster = rep(names(go)[i],nrow(tbl));
          cols = c("cluster.name","module.size","cluster.resolution","compactness.resolution","parent.cluster","compactness.parent","compactness.module",
                   "compactness.pvalue","intra.link.ratio","intra.link.pvalue")
          
          # rename modules according to the maximal module index number
          nm = max(as.integer(gsub("^M","",names(modules))))
          new.names = paste0("M",(nm+1):(nm + length(mods)));names(new.names) = names(mods)
          names(mods) = new.names[names(mods)]
          tbl$cluster.name = new.names[tbl$cluster.name]
          
          if (nrow(module.table) > 0)
          {
            mtbl = rbind.data.frame(module.table[,cols],tbl[,cols])
          }else{
            mtbl = rbind.data.frame(tbl[,cols])
          }
          
          comp.o = subset(mtbl,parent.cluster == "M0")$compactness.parent[1]
          gh = graph.data.frame(mtbl[,c("cluster.name","parent.cluster")],directed = TRUE)
          prnts = lapply(mtbl$cluster.name,function(x) setdiff(names(ego(graph = gh,order = 5,nodes = x,mode = "out")[[1]]),x))
          prnts.comp = sapply(prnts,function(x) {
            out = c(comp.o,subset(mtbl,cluster.name %in% x)$compactness.module)
            val = min(out,na.rm = T)
            return(val)
          })
          if (any(is.infinite(prnts.comp))) prnts.comp[is.infinite(prnts.comp)] = subset(mtbl,parent.cluster == "M0")$compactness.parent
          mtbl$compactness.parent.infimum = prnts.comp
          mtbl$is.compact = mtbl$compactness.module < mtbl$compactness.parent.infimum & mtbl$compactness.pvalue < pcut
          mtbl$is.coherent = mtbl$intra.link.pvalue < pcut
          if (valid.gate == "density") mtbl$is.valid = (mtbl$is.coherent) & mtbl$module.size >= min.size
          if (valid.gate == "compact") mtbl$is.valid = (mtbl$is.compact) & mtbl$module.size >= min.size
          if (valid.gate == "both") mtbl$is.valid = (mtbl$is.compact & mtbl$is.coherent) & mtbl$module.size >= min.size
          mtbl = subset(mtbl,cluster.name %in% names(mods))
          if (any(mtbl$is.valid))
          {
            # update modules
            modules = c(modules,mods)
            module.table = rbind.data.frame(module.table,mtbl)
            
            gc = lapply(mods[subset(mtbl,is.valid)$cluster.name],function(x) subgraph(go[[i]],vids = x))
            gn = c(gn,gc)
          }else{
            indiv.run[i] = FALSE
          }
        }else{
          indiv.run[i] = FALSE
        }
        
      }
  
      gn = gn[sapply(gn,vcount) >= min.size]
      if (any(indiv.run) & length(gn) > 0)
      {
        do.run = TRUE
        go = gn
      }else{
        do.run = FALSE
      }
      
    }
    
    output = list(modules = modules,module.table = module.table,cell.network = g.in,alpha.value = alpha)
    return(output)
  }

#' multi-scale clustering workflow
#' 
#' @description Wrapper function to perform LEN and top-down iterative clustering altogether. 
#' @docType package
#' @param dat A log-normalized single-cell expression matrix (rows: genes, columns: cell).
#' @param bcnt A binary matrix of genes (row) by cells (column). If expressed, \code{bcnt[i,j] = 1}, otherwise 0.   
#' @param n.cores An integer for the number of cores for parallelization 
#' @param min.cells An integer specifying the minumum cell cluster size
#' @param sim.func A function to calculate cell-cell similarity. Default is the pairwise Pearson's correlation 
#' @param is.decreasing A logical. TRUE if dist.func is similarity function to show greater value corresponds to greater similarity. FALSE if dist.func is dissimilarity.  
#' @param d.func A function to convert the pairwise similarity to pairwise distance.  
#' @param valid.gate Filtering method to identify significant clusters. "density" compares intra-cluster link density between parent and child clusters. "compact" compared compactness with the parent clusters. "both" uses "density" and "compact" filters simultaneously. Default is "both".
#' @return Returns a list containing clustering results with optimal.modules (the final results from adapt split), resolution.tested (clustering resolution sweeped), resolution.used (the final value of clustering resolution), and module.table containing the compactness and intra-cluster connectivity stats.
#' @export
#' @examples 
#'\donttest{
#' data(simMix1)
#' library(scater)
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
#' # for imputed data for correlation computation 
#' dat = subset(qc.out$gene.data.alra,bio > 0 & p.value < 0.05)
#' bcnt.alra = counts(qc.out$sce)[rownames(dat),]
#' m.alra = as.matrix(assay(qc.out$sce,"ALRA.logcounts")[rownames(dat),])
#' bcnt.alra[bcnt.alra > 1E-320] = 1
#' msc.res = msc_workflow(dat = m.alra,bcnt = bcnt.alra,n.cores = 4)
#' }
msc_workflow <- function(dat,bcnt,n.cores = 4,min.cells = 5,
                                  sim.func = function(a, b) cor(a,b,method ="pearson"),is.decreasing = TRUE,
                                  d.func = function(x) sqrt(2*(1-x)),
                                  valid.gate = "both")
  {
    ### get cell network
    gf = generate_cell_network.wt.loess(mat = dat,bcnt = bcnt,
                                        dist.func = sim.func,
                                        is.decreasing = is.decreasing,
                                        n.cores = n.cores)
    
    ## Run iterative clustering
    # identify scale parameter
    if (is.connected(gf))
    {
      alpha.res = identify_scale_parameter(g = gf,nv = 100,nr.max = 3,seed = 1234,d.func = d.func,mode = "diameter")
      gc()
      
      iter.res = iterative_clustering.par(g.in = gf,min.size = min.cells,d.func = d.func,alpha = alpha.res,
                                          n.cores = n.cores,valid.gate = valid.gate)
      gc()
      iter.res$cell.network = gf;
    }else{
      cout = components(gf)
      ci = which(cout$csize > log(vcount(gf)))
      glst = lapply(ci,function(n) induced_subgraph(graph = gf,vids = names(cout$membership)[cout$membership == n]))
      
      alpha.res = mean(sapply(glst,function(gi) identify_scale_parameter(g = gi,nv = 20,nr.max = 3,seed = 1234,d.func = d.func,mode = "diameter")))
      gc()
      
      iter.res = iterative_clustering.par(g.in = gf,min.size = min.cells,d.func = d.func,alpha = alpha.res,
                                          n.cores = n.cores,valid.gate = valid.gate)
    }
    iter.res$cell.network = gf;
    return(iter.res)
  }

#' Pruning cluster hierarchy
#' 
#' @description Quick filtering of module hierarchy tree to remove branches that do not anchor on the grand parent cluster, M0. 
#' @docType package
#' @param iter.res An output from iterative_clustering or iterative_clustering.par.
#' @return Adds "prunted.table" into the input data. 
#' @seealso [iterative_clustering()] or [iterative_clustering.par()] for the input. 
#' @name prune_hierarchy.compact
#' @export
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
