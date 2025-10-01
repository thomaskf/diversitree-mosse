## Models should provide:
##   1. make
##   2. info
##   3. make.cache, including initial tip conditions
##   4. initial.conditions(init, pars,t, idx)
##   5. rootfunc(res, pars, ...)

## Common other functions include:
##   stationary.freq
##   starting.point
##   branches

## 1: make
make.protracted.sampling <- function (tree,control=list(compiled=FALSE,backend="deSolve")) {
	cache <- make.cache.protracted.sampling(tree,control)
	all.branches <- make.all.branches.protracted.sampling(cache,cache$control)
	rootfunc <- function (res, pars, condition.surv=TRUE, root=ROOT.FLAT, root.p=NULL) {
  		vals <- res$vals
  		idx <- res$idx
  		n <- pars[4]
  		lq <- res$lq
  		d.root <- vals[c(4,(7+n+n):length(vals))]
  		if ( condition.surv ) { #condition the likelihood on that the survivial of the root
    		b <- pars[1]
    		#if (n>0) {
    		#	e.root <- vals[c(2,c(1,2,7:(6+n))[idx])]
    		#} else {
    		#	e.root <- vals[2]
    		#}
    		d.root[1] <- d.root[1] / (b * (1 - vals[2])^2)
    		d.root[-1] <- d.root[-1] / (b * (1 - vals[2]) * (1 - vals[1]))
  		}
  		root.p <- root.p.calc(d.root, pars, root, root.p, root.equi=F)
    	loglik <- log(sum(root.p * d.root)) + sum(lq)

  		loglik
  	}
	ll <- function (pars, condition.surv=TRUE) {
		pars <- c(pars,length(which(tree$sampling.f<1)))
		ans <- all.branches(pars)
		rootfunc(ans, pars, condition.surv)
	}
	class(ll) <- c("protracted","dtlik","function")
	ll
}

## 3: make.cache (& initial conditions)
make.cache.protracted.sampling <- function (tree,control) {
	## 1: tree
	cache <- make.cache.tree.protracted(tree)
  	cache$sampling.f <- tree$sampling.f
  	
  	## 2: Control structure
	control <- check.control.ode(control)
	cache$control <- control
	cache$control$h <- 1
	
	#info
	cache$info <- make.info.protracted.sampling(tree)
	
	cache
}

initial.tip.protracted.sampling <- function (cache) {
	group <- cache$group # a list with eeach cell corresponding to a clade, which records the order of branches to calculate within the clade
	n.tip <- cache$n.tip # total number of tips in the tree
	len <- cache$len # branch length of each branch
	sampling.f <- cache$sampling.f
	tmp0 <- which(sampling.f==1)
	tmp1 <- which(is.na(sampling.f))
	tmp2 <- which(sampling.f<1)
	sampling.f <- sampling.f[tmp2]
	#init is a function that generates y=all possible scenarios of the tip initial states in a clade, tip=idx of tip branches, t=tip branch length
	init <- function (group2,n.tip,len,sampling.f,tmp0,tmp1,tmp2) {
		tip <- group2[group2<=n.tip] # find the tip branches
		# the initial state is y, which is a matrix with rows are tips, columns are c(EI,EG,DI,DR)
		if (length(tip)>0) {
		t <- len[tip]
		k <- length(tip)
		if (k==1) {
			if (is.element(tip,tmp0)) {
				idx <- 1
				Drsp <- 1
			} else if (is.element(tip,tmp1)) {
				idx <- 2
				Drsp <- 1
			} else {
				idx <- which(tmp2==tip)
				Drsp <- sampling.f[idx]
				idx <- idx+2
			}
			y <- c(0,1,0,0,0,0,1-sampling.f,numeric(length(sampling.f)),Drsp)
		} else {
			if (is.element(tip[1],tmp0)) {
				idx <- 1
				Drsp <- 1
				y <- vector("list",length(tip)) 
				for (i in 2:length(tip)) {
					y[[i]] <- c(0,1,0,0,1,0,1-sampling.f,numeric(length(sampling.f)),0)
				}
			} else if (is.element(tip[1],tmp1)) {
				idx <- 2
				Drsp <- 1
				y <- vector("list",length(tip)) 
				for (i in 2:length(tip)) {
					y[[i]] <- c(0,1,0,0,0,1,1-sampling.f,numeric(length(sampling.f)),0)
				}
			} else {
				idx <- which(tmp2==tip[1])
				Drsp <- sampling.f[idx]
				tmp <- numeric(length(sampling.f))
				tmp[idx] <- 1
				idx <- idx+2
				y <- vector("list",length(tip)) 
				for (i in 2:length(tip)) {
					y[[i]] <- c(0,1,0,0,0,0,1-sampling.f,tmp,0)
				}
			}
			y[[1]] <- c(0,1,0,0,0,0,1-sampling.f,numeric(length(sampling.f)),Drsp)
		}
			list(y=y,tip=tip,t=t,idx=idx)
		} else {
			list(NULL)
		}
	}
	lapply(group,init,n.tip,len,sampling.f,tmp0,tmp1,tmp2) #apply init to each cell in the group list.
}

initial.conditions.protracted.sampling <- function(init, pars) {
	n <- pars[4]
	idx1 <- init[[1]]$idx
	idx2 <- init[[2]]$idx
	y1 <- init[[1]]$y
	y2 <- init[[2]]$y
	tmp <- which(y1[-c(1:(6+n+n))]==0)
	if (length(tmp)>0) {
		y1 <- y1[c(1:(6+n+n),c((7+n+n):length(y1))[-tmp])]
		idx1 <- idx1[-tmp]
	}
	tmp <- which(y2[-c(1:(6+n+n))]==0)
	if (length(tmp)>0) {
		y2 <- y2[c(1:(6+n+n),c((7+n+n):length(y1))[-tmp])]
		idx2 <- idx2[-tmp]
	}
	if (n>0) {
		y <- c(y1[1:2],
	  		y1[3] * y2[3] * pars[1],
	  		y1[4] * y2[3] * pars[1] + y1[3] * y2[4] * pars[1],
	  		y1[5] * y2[5] * pars[1],
	  		y1[6] * y2[6] * pars[1],
	  		y1[7:(6+n)],
	  		y1[(7+n):(6+n+n)] * y2[(7+n):(6+n+n)] * pars[1],
	  		y1[-c(1:(6+n+n))] * c(y2[5:6],y2[(7+n):(6+n+n)])[idx1] * pars[1],
	  		c(y1[5:6],y1[(7+n):(6+n+n)])[idx2] * y2[-c(1:(6+n+n))] * pars[1])
	} else {
		y <- c(y1[1:2],
	  		y1[3] * y2[3] * pars[1],
	  		y1[4] * y2[3] * pars[1] + y1[3] * y2[4] * pars[1],
	  		y1[5] * y2[5] * pars[1],
	  		y1[6] * y2[6] * pars[1],
	  		y1[-c(1:(6+n+n))] * c(y2[5:6],y2[(7+n):(6+n+n)])[idx1] * pars[1],
	  		c(y1[5:6],y1[(7+n):(6+n+n)])[idx2] * y2[-c(1:(6+n+n))] * pars[1])
	}
	idx <- c(idx1,idx2)
	list(y=y,idx=idx)
}

######################################################################
## Extra core stuff:
#this is similar to make.all.branches.dtlik
make.all.branches.protracted.sampling <- function (cache,control) {
	branches <- make.branches.protracted.sampling(cache$info,control)
	function (pars) {
		cache$y <- initial.tip.protracted.sampling(cache)
		all.branches.matrix.protracted.sampling(pars,cache,initial.conditions=initial.conditions.protracted.sampling,branches)
	}
}

make.branches.protracted.sampling <- function(info,control) {
	eps <- control$eps
	ode <- list(make.ode.sampling(info[[1]],control),make.ode.sampling(info[[2]],control)) #make a list of ODEs, with the first cell for DI* and the second cell for DI.
	branches1 <- function(y, len, pars, t0, idx, star) {ode[[star]](y, c(t0, t0+len), c(pars,idx))} #function to integrate ODE along a branch
  	make.branches.comp.sampling(branches1, eps)
}

make.ode.sampling <- function(info, control) {
  #control <- check.control.ode(control)
  #info <- check.info.ode(info, control)
  backend  <- control$backend
  if ( backend == "gslode" )
    ode <- make.ode.gslode(info, control)
  else if ( backend == "deSolve" )
    ode <- make.ode.deSolve(info, control)
  else # should have been prevented by now
    stop("Invalid backend", backend)
  ode
}

all.branches.matrix.protracted.sampling <- function (pars,cache,initial.conditions,branches) {
	len <- cache$len #branch length
  	depth <- cache$depth #branch starting time
  	children <- cache$children #children of each node
  	order <- cache$order[-length(cache$order)] #parent of each node
  	root <- cache$root #root label
	group <- cache$group #list of clades
	stars <- cache$stars #indicator for which ODE to use
	h <- cache$control$h
	
  	n <- length(len)
  	lq <- rep(0, n)
 	n.tip <- cache$n.tip

  	y <- cache$y #initial tip states
  	branch.init <- branch.base <- vector("list", n)
	
	#start from each clade
	for ( x in 1:length(y) ) {
        idx <- y[[x]]$tip #tip index in each clade
        if (length(idx)==1) { #for clade with single tip
        	branch.init[[idx]] <- y[[x]]$y
        	ans <- branches(y[[x]]$y, y[[x]]$t, pars, 0, y[[x]]$idx, star=2)
        	lq[idx] <- ans[[1]]
        	branch.base[[idx]] <- list(y=ans[[2]],idx=y[[x]]$idx)
        } else if (length(idx)>1){ #for clade more than one tip
        	idx2 <- group[[x]]
        	idx2 <- idx2[idx2>n.tip]
        	#calculate along tip branch
        	for (j in 1:length(idx)) { 
        		branch.init[[idx[j]]] <- y[[x]]$y[[j]]
        		ans <- branches(y[[x]]$y[[j]], y[[x]]$t[j], pars, 0, y[[x]]$idx, star=1) #use the ODE for DI*, so star=1
         		lq[idx[j]] <- ans[[1]]
        		branch.base[[idx[j]]] <- list(y=ans[[2]],idx=y[[x]]$idx)
        	}
        	#calculate along internal branch
         	for (j in idx2[-length(idx2)]) {
         		tmp <- children[j,!is.element(children[j,],group[[x]])] 
         		if (length(tmp)>0) {
         			if (is.matrix(branch.base[[tmp]])) {
         				branch.base[[tmp]][,(h+2):(2*h+1)] <- 0 ###need to update
         			} else {
         				branch.base[[tmp]][h+3] <- 0
         				branch.base[[tmp]][-c(1:(6+n+n))] <- 0
         			}
         		} #find children that leads to another species and set its DR=0
        		y.in <- initial.conditions(branch.base[children[j,]], pars) #calculate the boundary condition of the branch
        		if ( !all(is.finite(y.in$y)) )
      				stop("Bad initial conditions: calculation failure along branches?")
    			branch.init[[j]] <- y.in$y
    			ans <- branches(y.in$y, len[j], pars, depth[j], y.in$idx, star=1)
    			lq[j] <- ans[[1]]
    			branch.base[[j]] <- list(y=ans[[2]],idx=y.in$idx)
        	}
        	#calculate the boundary conditions for the basal branch
        	j <- idx2[length(idx2)]
        	tmp <- children[j,!is.element(children[j,],group[[x]])] 
         	if (length(tmp)>0) {
         		if (is.matrix(branch.base[[tmp]])) {
         				branch.base[[tmp]][,(h+2):(2*h+1)] <- 0
         			} else {
         				branch.base[[tmp]][h+3] <- 0
         				branch.base[[tmp]][-c(1:(6+n+n))] <- 0
         			}
         	} #find children that leads to another species and set its DR=0
        	y.in <- initial.conditions(branch.base[children[j,]], pars)
        	if ( !all(is.finite(y.in$y)) )
      			stop("Bad initial conditions: calculation failure along branches?")
    		branch.init[[j]] <- y.in$y
    		#calculate along the basal branch
    		if (j!=root) {
        		ans <- branches(y.in$y, len[j], pars, depth[j], y.in$idx, star=2) #use the ODE for DI, so star=2
        		lq[j] <- ans[[1]]
        		branch.base[[j]] <- list(y=ans[[2]],idx=y.in$idx)
        	}
        } else { #if the group is not a clade, but the branch linking a clade to another clade that it is nested in, then we need to calculate this branch first.
        	for (j in group[[x]]) {
        		y.in <- initial.conditions(branch.base[children[j,]], pars)
        		if ( !all(is.finite(y.in$y)) )
      				stop("Bad initial conditions: calculation failure along branches?")
    			branch.init[[j]] <- y.in$y
    			if (j!=root) {
    				ans <- branches(y.in$y, len[j], pars, depth[j], y.in$idx, star=2)
    				lq[j] <- ans[[1]]
        			branch.base[[j]] <- list(y=ans[[2]],idx=y.in$idx)
        		}
        	}
        }
    }
	for ( i in order ) {
		y.in <- initial.conditions(branch.base[children[i,]], pars)
        if ( !all(is.finite(y.in$y)) )
      		stop("Bad initial conditions: calculation failure along branches?")
      	branch.init[[i]] <- y.in$y
      	if (i!=root) {
    		ans <- branches(y.in$y, len[i], pars, depth[i], y.in$idx, star=2)
    		lq[i] <- ans[[1]]
    		branch.base[[i]] <- list(y=ans[[2]],idx=y.in$idx)
		}
  	}

	y.in <- initial.conditions(branch.base[children[root,]], pars)
   	branch.init[[root]] <- y.in$y
	list(init=branch.init, base=branch.base, lq=lq, vals=y.in$y, idx=y.in$idx)
}   	

make.branches.comp.sampling<- function(branches1, eps=0) {
    function(y, len, pars, t0, idx, star) {
      ret <- branches1(y, len, pars, t0, idx, star)
      n <- pars[4]
      #comp.idx <- c(4,(7+n+n):length(y))
      #q <- sum(ret[comp.idx])
      q <- 10^-3
      if (n>0) {
      	comp.idx <- c(3:6,(7+n):length(y))
      } else {
      	comp.idx <- c(3:6)
      }
      if ( all(ret[comp.idx] >= eps) && all(ret[comp.idx] <= 1) ) {
      	if ( any(ret[comp.idx] <= 10^-5 && ret[comp.idx] >0) ) {
        	ret[comp.idx] <- ret[comp.idx] / q
        	lq <- log(q)
        } else {
        	lq <- 0
        }
        list(lq, ret)
      } else {
        ti <- len[length(len)]/2
        len1 <- c(len[len <= ti], ti)
        len2 <- len[len > ti] - ti
        n1 <- length(len1)

        ret1 <- Recall(y, len1, pars, t0, idx, star)
        ret2 <- Recall(ret1[[2]][,n1], len2, pars, t0 + ti, idx, star)
        ret2[[1]] <- ret2[[1]] + ret1[[1]][n1]

        list(c(ret1[[1]][-n1], ret2[[1]]),
             cbind(ret1[[2]][,-n1], ret2[[2]]))
      }
    }
}

make.info.protracted.sampling <- function (phy) {
	list(
	list(name="protracted.sampling",
	     derivs=derivs.protracted.sampling.star,
		 np=3L,
		 ny=NA,
		 k=NA,
		 idx.e=NA,
		 idx.d=NA,
		 argnames=default.argnames.protracted(),
		 phy=phy,
		 ml.default="subplex",
		 mcmc.lowerzero=TRUE),
	list(name="protracted.sampling",
	     derivs=derivs.protracted.sampling,
		 np=3L,
		 ny=NA,
		 k=NA,
		 idx.e=NA,
		 idx.d=NA,
		 argnames=default.argnames.protracted(),
		 phy=phy,
		 ml.default="subplex",
		 mcmc.lowerzero=TRUE)
	)	 
}

default.argnames.protracted <- function() {
	c("b","mu","lambda")
}

#E(Esp0),Er,Di,Dr,Disp0,Disp1,Esp,Disp,Drsp*
derivs.protracted.sampling.star <- function (t,y,pars) {
			n <- pars[4]
			Ei <- y[1]
			Er <- y[2]
			Di <- y[3]
			Dr <- y[4]
			Disp0 <- y[5]
			Disp1 <- y[6]
			if (n>0) {
				Esp <- y[7:(6+n)]
				Disp <- y[(7+n):(6+n+n)]
				E <- c(Ei,Er,Esp)
			} else {
				E <- c(Ei,Er)
			}
			Drsp <- y[-c(1:(6+n+n))]
			tmp <- sum(Drsp)
			b <- pars[1]
			mu <- pars[2]
			lambda <- pars[3]
			idx <- pars[-c(1:4)]
			N <- length(idx)
			res <- y
			res[1:6] <- c(-(b+mu)*Ei + b*Ei*Ei + mu,
			  -(b+mu+lambda)*Er + b*Er*Er + mu + lambda*Ei,
			  0,
			  0,
			  -(b+mu+lambda)*Disp0 + 2*b*Disp0*Ei,
			  -(b+mu+lambda)*Disp1 + 2*b*Disp1*Er)
			if (n>0) {
				res[7:(6+n)] <- -(b+mu)*Esp + b*Esp*Esp + mu
				res[(7+n):(6+n+n)] <- -(b+mu+lambda)*Disp + 2*b*Disp*Esp
			}
			res[-c(1:(6+n+n))] <- -(b+mu+lambda)*Drsp + 2*b*Drsp*E[idx]
			res
}

derivs.protracted.sampling <- function (t,y,pars) {
			n <- pars[4]
			Ei <- y[1]
			Er <- y[2]
			Di <- y[3]
			Dr <- y[4]
			Disp0 <- y[5]
			Disp1 <- y[6]
			if (n>0) {
				Esp <- y[7:(6+n)]
				Disp <- y[(7+n):(6+n+n)]
				E <- c(Ei,Er,Esp)
			} else {
				E <- c(Ei,Er)
			}
			Drsp <- y[-c(1:(6+n+n))]
			tmp <- sum(Drsp)
			b <- pars[1]
			mu <- pars[2]
			lambda <- pars[3]
			idx <- pars[-c(1:4)]
			N <- length(idx)
			res <- y
			res[1:6] <- c(-(b+mu)*Ei + b*Ei*Ei + mu,
			  -(b+mu+lambda)*Er + b*Er*Er + mu + lambda*Ei,
			  -(b+mu+lambda)*Di + 2*b*Di*Ei + lambda*Dr + lambda*tmp,
			  -(b+mu)*Dr + 2*b*Dr*Ei + lambda*tmp,
			  -(b+mu+lambda)*Disp0 + 2*b*Disp0*Ei + lambda*Dr + lambda*tmp,
			  -(b+mu+lambda)*Disp1 + 2*b*Disp1*Er + lambda*Dr + lambda*tmp)
			if (n>0) {
				res[7:(6+n)] <- -(b+mu)*Esp + b*Esp*Esp + mu
				res[(7+n):(6+n+n)] <- -(b+mu+lambda)*Disp + 2*b*Disp*Esp + lambda*Dr + lambda*tmp
			}
			res[-c(1:(6+n+n))] <- -(b+mu+lambda)*Drsp + 2*b*Drsp*E[idx]
			res
}