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
make.protracted.hybrid <- function (tree,control=list(compiled=FALSE,backend="deSolve")) {
	cache <- make.cache.protracted.hybrid(tree,control)
	all.branches <- make.all.branches.protracted.hybrid(cache,cache$control)
	rootfunc <- function (res, pars, t=max(cache$depth), condition.surv) {
  		vals <- res$vals
  		lq <- res$lq
  		d.root <- vals[3] #the root must be in state R
  
  		if ( condition.surv ) { #condition the likelihood on that the survivial of the root
  			#t <- tree$root.depth
    		b <- pars[1]
    		e.root <- vals[1]
    		func <- function (t,y,pars) {
    			b <- pars[1]
    			mu <- pars[2]
    			lambda <- pars[3]
    			Ei <- min(y[1],1)
    			Er <- min(y[2],1)
    			list(c(-(b+mu)*Ei + b*Ei*Ei + mu,
    			-(b+mu+lambda)*Er + b*Er*Er + mu + lambda * Ei))
    		}
    		er.root <- ode(y=c(0,1),times=c(0,t),func=func,pars)[2,2]
    		d.root <- d.root / (b * (1 - e.root) * (1 - er.root))
  		}

  		loglik <- log(d.root) + sum(lq)

  		loglik
  	}
	ll <- function (pars, condition.surv=T) {
		ans <- all.branches(pars)
		rootfunc(ans, pars, max(cache$depth), condition.surv)
	}
	class(ll) <- c("protracted","dtlik","function")
	ll
}

make.cache.protracted.hybrid <- function (tree,control) {
	## 1: tree
  	#tree <- check.tree(tree)
  	
  	## 2: Control structure
	cache <- make.cache.tree.protracted(tree)
	control <- check.control.ode(control)
	cache$control <- control
	cache$control$h <- 1
	
	cache$info <- make.info.protracted.representative(tree)
	
	cache
}

initial.conditions.protracted.hybrid <- function(init, pars, t)
  c(init[[1]][1],init[[1]][2],init[[1]][3] * init[[2]][3] * pars[1],
    (init[[1]][3] * init[[2]][4] + init[[1]][4] * init[[2]][3]) * pars[1])

initial.tip.protracted.hybrid <- function (cache) {
	group <- cache$group # a list with eeach cell corresponding to a clade, which records the order of branches to calculate within the clade
	n.tip <- cache$n.tip # total number of tips in the tree
	len <- cache$len # branch length of each branch
	#init is a function that generates y=all possible scenarios of the tip initial states in a clade, target=idx of tip branches, t=tip branch length
	init <- function (group2,n.tip,len) {
		target <- group2[group2<=n.tip] # find the tip branches
		# the initial state is y, which is a matrix with rows are tips, columns are c(EI,EG,DI,DR)
		if (length(target)>0) {
		t <- len[target]
		k <- length(target)
		if (k==1) {
			y <- c(0,1,0,1)
		} else {
			y <- rep(list(c(0,1,1,0)),k)
			y[[1]] <- c(0,1,0,1)
		}
			list(y=y,target=target,t=t)
		} else {
			list(NULL)
		}
	}
	lapply(group,init,n.tip,len) #apply init to each cell in the group list.
}

######################################################################
## Extra core stuff:
#this is similar to make.all.branches.dtlik
make.all.branches.protracted.hybrid <- function (cache,control) {
	branches1 <- make.branches.protracted.representative(cache$info,cache$control)
	branches2 <- function(y, len, pars, t0, star) {
		b <- pars[1]
		mu <- pars[2]
		lambda <- pars[3]
		r <- b-mu
		z <- exp(len * r)
		z2 <- exp(- len * lambda)
		e0 <- y[1]
  		dr0 <- y[4]
  		y[1] <- (mu + z*(e0 - 1)*mu - b*e0) / (mu + z*(e0 - 1)*b - b*e0)
  		A <- (z * r * r)/(z * b - mu + (1-z)*b*e0)^2
  		if (star==1) {
  			y[-c(1:2)] <- A * y[-c(1:2)] * z2
  		} else if (star==2) {
  			y[4] <- A * y[4]
  			y[3] <- y[4] + (y[3]-dr0) * A * z2
  		}
  		if (any(y[-c(1:2)]<=10^-5 && y[-c(1:2)]>0)) {
  			lq <- sum(y[-c(1:2)])
  			y[-c(1:2)] <- y[-c(1:2)] / lq
  			lq <- log(lq)
  		} else {
  			lq <- 0
  		}
  		list(lq=lq, y=y) 		
		}

	function (pars) {
		cache$y <- initial.tip.protracted.hybrid(cache)
		all.branches.matrix.protracted.hybrid(pars,cache,initial.conditions=initial.conditions.protracted.hybrid,branches1,branches2)
	}
}

all.branches.matrix.protracted.hybrid <- function (pars,cache,initial.conditions,branches1,branches2) {
	len <- cache$len #branch length
  	depth <- cache$depth #branch starting time
  	children <- cache$children #children of each node
  	order <- cache$order #parent of each node
  	root <- cache$root #root label
	group <- cache$group #list of clades
	stars <- cache$stars #indicator for which ODE to use
	h <- cache$control$h
	
  	n <- length(len)
  	lq <- rep(0, n)  #this is the log of DI+DR after the calculation on each branch. It is used to standardize DI and DR to avoid rounding error due to extremely small likelihood value
 	n.tip <- cache$n.tip

  	y <- cache$y #initial tip states
  	branch.init <- branch.base <- vector("list", n)
	
	#start from each clade
    for ( x in 1:length(y) ) {
        idx <- y[[x]]$target #tip index in each clade
        if (length(idx)==1) { #for clade with single tip
        	branch.init[idx] <- list(y[[x]]$y)
        	ans <- branches1(y[[x]]$y, y[[x]]$t, pars, 0)
        	lq[idx] <- ans[[1]]
        	branch.base[[idx]] <- ans[[2]]
        } else if (length(idx)>1){ #for clade more than one tip
        	idx2 <- group[[x]]
        	idx2 <- idx2[idx2>n.tip]
        	#calculate along tip branch
        	for (j in 1:length(idx)) { 
        		branch.init[idx[j]] <- y[[x]]$y[j]
        		ans <- branches2(y[[x]]$y[[j]], y[[x]]$t[j], pars, 0, star=1) #use the ODE for DI*, so star=1
         		lq[idx[j]] <- ans[[1]]
        		branch.base[[idx[j]]] <- ans[[2]]
        	}
        	#calculate along internal branch
         	for (j in idx2[-length(idx2)]) {
         		tmp <- children[j,!is.element(children[j,],group[[x]])] 
         		if (length(tmp)>0) {
         			if (is.matrix(branch.base[[tmp]])) {
         				branch.base[[tmp]][,(h+2):(2*h+1)] <- 0
         			} else {
         				branch.base[[tmp]][h+2] <- 0
         			}
         		} #find children that leads to another species and set its DR=0
        		y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j]) #calculate the boundary condition of the branch
        		if ( !all(is.finite(y.in)) )
      				stop("Bad initial conditions: calculation failure along branches?")
    			branch.init[[j]] <- y.in
    			ans <- branches2(y.in, len[j], pars, depth[j], star=1)
    			lq[j] <- ans[[1]]
    			branch.base[[j]] <- ans[[2]]
        	}
        	#calculate the boundary conditions for the basal branch
        	j <- idx2[length(idx2)]
        	tmp <- children[j,!is.element(children[j,],group[[x]])] 
         	if (length(tmp)>0) {
         		if (is.matrix(branch.base[[tmp]])) {
         				branch.base[[tmp]][,(h+2):(2*h+1)] <- 0
         			} else {
         				branch.base[[tmp]][h+2] <- 0
         			}
         	} #find children that leads to another species and set its DR=0
        	y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j])
        	if ( !all(is.finite(y.in)) )
      			stop("Bad initial conditions: calculation failure along branches?")
    		branch.init[[j]] <- y.in
    		#calculate along the basal branch
    		if (j!=root) {
        		ans <- branches2(y.in, len[j], pars, depth[j], star=2) #use the ODE for DI, so star=2
        		lq[j] <- ans[[1]]
        		branch.base[[j]] <- ans[[2]]
        	}
        } else { #if the group is not a clade, but the branch linking a clade to another clade that it is nested in, then we need to calculate this branch first.
        	for (j in group[[x]]) {
        		y.in <- initial.conditions(branch.base[children[j,]], pars, depth[j])
        		if ( !all(is.finite(y.in)) )
      				stop("Bad initial conditions: calculation failure along branches?")
    			branch.init[[j]] <- y.in
    			if (j!=root) {
    				ans <- branches2(y.in, len[j], pars, depth[j], star=2)
    				lq[j] <- ans[[1]]
        			branch.base[[j]] <- ans[[2]]
        		}
        	}
        }
    }

    for ( i in order ) {
    	if (i!=root) {
    	y.in <- initial.conditions(branch.base[children[i,]], pars, depth[i])
    	if ( !all(is.finite(y.in)) )
      		stop("Bad initial conditions: calculation failure along branches?")
    	branch.init[[i]] <- y.in
    	ans <- branches2(y.in, len[i], pars, depth[i], star=2)
    	lq[i] <- ans[[1]]
    	branch.base[[i]] <- ans[[2]]
    	}
  	}
 	
  	y.in <- initial.conditions(branch.base[children[root,]], pars, depth[root])
  	branch.init[[root]] <- y.in
  	ans <- branches2(y.in, len[root], pars, depth[root], star=2)
    lq[root] <- ans[[1]]
    branch.base[[root]] <- ans[[2]]
  	
  list(init=branch.init, base=branch.base, lq=lq, vals=branch.base[[root]])
}

