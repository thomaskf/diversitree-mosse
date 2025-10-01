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
	rootfunc <- function (res, pars, condition.surv=TRUE, root=ROOT.OBS, root.p=NULL) {
  		vals <- res$vals
  		idx <- res$idx
  		n <- pars[4]
  		lq <- res$lq
  		d.root <- vals[(7+n+n):length(vals)]
  		d.root <- c(vals[4]-sum(d.root),d.root)
  		if ( condition.surv ) { #condition the likelihood on that the survivial of the root
    		b <- pars[1]
    		if (n>0) {
    			e.root <- vals[c(2,c(1,2,7:(6+n))[idx])]
    		} else {
    			e.root <- vals[2]
    		}
    		d.root <- d.root / (b * (1 - e.root[1]) * (1 - e.root))
  		}
  		root.p <- root.p.calc(d.root, pars, root, root.p, root.equi=F)

  		if ( root == ROOT.ALL )
    		loglik <- log(d.root) + sum(lq)
  		else
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
  	#tree <- check.tree(tree)
  	
  	## 2: Control structure
	cache <- make.cache(tree)
	control <- check.control.ode(control)
	cache$control <- control
	
	cache$info <- make.info.protracted.sampling(tree)
	cache$sampling.f <- tree$sampling.f
	
	cache
}

initial.tip.protracted.sampling <- function (cache) {
	tips <- cache$tips
	sampling.f <- cache$sampling.f
	tmp0 <- which(sampling.f==1)
	tmp1 <- which(is.na(sampling.f))
	tmp2 <- which(sampling.f<1)
	sampling.f <- sampling.f[tmp2]
	init <- function (tip,sampling.f,tmp0,tmp1,tmp2) {
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
		list(y=c(0,1,0,Drsp,0,0,1-sampling.f,numeric(length(sampling.f)),Drsp),idx=idx)
	}
	lapply(tips,init,sampling.f,tmp0,tmp1,tmp2) 
}

initial.conditions.protracted.sampling <- function(init, pars) {
	n <- pars[4]
	idx1 <- init[[1]]$idx
	idx2 <- init[[2]]$idx
	y1 <- init[[1]]$y
	y2 <- init[[2]]$y
	y1[4] <- y1[4]-sum(y1[-c(1:(6+n+n))])
	y2[4] <- y2[4]-sum(y2[-c(1:(6+n+n))])
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
	y[4] <- y[4]+sum(y[-c(1:(6+n+n))])
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
	ode <- make.ode.sampling(info,control) #make a list of ODEs, with the first cell for DI* and the second cell for DI.
	branches1 <- function(y, len, pars, t0, idx) {ode(y, c(t0, t0+len), c(pars,idx))} #function to integrate ODE along a branch
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
	
  	n <- length(len)
  	lq <- rep(0, n)
 	n.tip <- cache$n.tip

  	y <- cache$y #initial tip states
  	branch.init <- branch.base <- vector("list", n)
	
	#start from each clade
    for ( x in 1:n.tip ) {
    	branch.init[[x]] <- y[[x]]
    	ans <- branches(branch.init[[x]]$y, len[x], pars, 0, y[[x]]$idx)
    	lq[x] <- ans[[1]]
    	branch.base[[x]] <- list(y=ans[[2]],idx=y[[x]]$idx)
    }
    
	for ( i in order ) {
		y.in <- initial.conditions(branch.base[children[i,]], pars)
        if ( !all(is.finite(y.in$y)) )
      		stop("Bad initial conditions: calculation failure along branches?")
      	branch.init[[i]] <- y.in
    	ans <- branches(branch.init[[i]]$y, len[i], pars, depth[i], branch.init[[i]]$idx)
    	lq[i] <- ans[[1]]
    	branch.base[[i]] <- list(y=ans[[2]],idx=branch.init[[i]]$idx)
  	}
	
	y.in <- initial.conditions(branch.base[children[root,]], pars)
   	branch.init[[root]] <- y.in
	list(init=branch.init, base=branch.base, lq=lq, vals=y.in$y, idx=y.in$idx)
}   	

make.branches.comp.sampling<- function(branches1, eps=0) {
    function(y, len, pars, t0, idx) {
      ret <- branches1(y, len, pars, t0, idx)
      n <- pars[4]
      comp.idx <- c(4,(7+n+n):length(y))
      q <- sum(ret[comp.idx])
      comp.idx <- c(3:6,(7+n):length(y))
      if ( all(ret[comp.idx] >= eps) ) {
        ret[comp.idx] <- ret[comp.idx] / q
        lq <- log(q)
        list(lq, ret)
      } else {
        ti <- len[length(len)]/2
        len1 <- c(len[len <= ti], ti)
        len2 <- len[len > ti] - ti
        n1 <- length(len1)

        ret1 <- Recall(y, len1, pars, t0, idx)
        ret2 <- Recall(ret1[[2]][,n1], len2, pars, t0 + ti, idx)
        ret2[[1]] <- ret2[[1]] + ret1[[1]][n1]

        list(c(ret1[[1]][-n1], ret2[[1]]),
             cbind(ret1[[2]][,-n1], ret2[[2]]))
      }
    }
}

make.info.protracted.sampling <- function (phy) {
	list(name="protracted.samoling",
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
}

#E(Esp0),Er,Di,Dr,Disp0,Disp1,Esp,Disp,Drsp*
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
			b <- pars[1]
			mu <- pars[2]
			lambda <- pars[3]
			idx <- pars[-c(1:4)]
			#N <- length(idx)
			res <- y
			res[1:6] <- c(-(b+mu)*Ei + b*Ei*Ei + mu,
			  -(b+mu+lambda)*Er + b*Er*Er + mu + lambda*Ei,
			  -(b+mu+lambda)*Di + 2*b*Di*Ei + lambda*Dr,
			  -(b+mu)*Dr + 2*b*Dr*Ei,
			  -(b+mu+lambda)*Disp0 + 2*b*Disp0*Ei + lambda*Dr,
			  -(b+mu+lambda)*Disp1 + 2*b*Disp1*Er + lambda*Dr)
			if (n>0) {
				res[7:(6+n)] <- -(b+mu)*Esp + b*Esp*Esp + mu
				res[(7+n):(6+n+n)] <- -(b+mu+lambda)*Disp + 2*b*Disp*Esp + lambda*Dr
			}
			res[-c(1:(6+n+n))] <- -(b+mu+lambda)*Drsp + 2*b*Drsp*E[idx]
			res
}