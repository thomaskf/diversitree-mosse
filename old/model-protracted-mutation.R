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
make.protracted.mutation <- function (tree, states, states.sd, b, lambda, mu, control=NULL) {
	#the tree has edge.spec is 0/1 recording wether the edge (index by its descendant node) has a speciation event. c(node.label, tip.label)
	cache <- make.cache.protracted.mutation(tree, states, states.sd, b, mu, lambda, control)
	all.branches <- make.all.branches.protracted.mutation(cache, cache$control)
	rootfunc <- make.rootfunc.protracted.mutation(cache)
	f.pars <- make.pars.protracted.mutation(cache)
	ll <- function (pars, condition.surv=TRUE, root=ROOT.OBS, root.f=NULL) {
		pars2 <- f.pars(pars)
		ans <- all.branches(pars2)
		rootfunc(ans, pars2, condition.surv, root, root.f)
	}
	class(ll) <- c("protracted","dtlik","function")
	ll
}

## 2: info
make.info.protracted.mutation <- function(b, mu, lambda, phy) {
  ## Work around for .split:
  if ( !is.null(lambda) )
    argnames <- default.argnames.protracted.mutation(b, mu, lambda)
  else
    argnames <- NULL
  list(name="protracted",
       ## Parameters:
       np=NA,
       argnames=argnames,
       ## Variables:
       ny=NA,
       k=NA,
       idx.e=NA,
       idx.d=NA,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=FALSE, # not for many models.
       ## These are optional
       doc=NULL,
       reference=NA,
       ## These are special to ProSSE:
       b=b,
       lambda=lambda,
       mu=mu)
}
default.argnames.protracted.mutation <- function(b, mu, lambda) {
  c(sprintf("b.%s", names(formals(b))[-1]),
    sprintf("m.%s", names(formals(mu))[-1]),
    sprintf("l.%s", names(formals(lambda))[-1]),
    "drift", "diffusion")
}

## 3: make.cache (& initial conditions)
make.cache.protracted.mutation <- function (tree, states, states.sd, b, mu, lambda, control) {
	## 1: tree
  	#tree <- check.tree(tree)
  	
  	## 2: states & errors
  	tmp <- check.states.quasse(tree, states, states.sd)
  	states <- tmp$states
 	states.sd <- tmp$states.sd
 	
 	## 3: Control structure (lots of checking!)
 	control <- check.control.quasse(control, tree, states)
 	
	cache <- make.cache.tree.protracted(tree)
	cache$states  <- states
    cache$states.sd <- states.sd
	cache$control <- control
	cache$control$h <- 1
	
	## 4: Speciation/extinction functions
    n.b <- check.f.quasse(b)
    n.mu     <- check.f.quasse(mu)
    n.lambda <- check.f.quasse(lambda)
    n.args   <- n.b + n.lambda + n.mu + 2
    args <- list(b=seq_len(n.b),
                 mu=seq_len(n.mu) + n.b,
                 lambda=seq_len(n.lambda) + n.b + n.mu,
                 drift=n.b + n.mu + n.lambda + 1,
                 diffusion=n.b + n.mu + n.lambda + 2)
	
	cache$b <- b
    cache$lambda <- lambda
    cache$mu <- mu
    cache$args <- args

	cache$info <- make.info.protracted.mutation(b, mu, lambda, tree)
	
	cache
}

## 4: initial.conditions
initial.tip.protracted.mutation <- function (cache, control, x) {
	nx <- control$nx * control$r
    npad <- nx - length(x)
    states <- cache$states
    states.sd <- cache$states.sd
	group <- cache$group # a list with each cell corresponding to a clade, which records the order of branches to calculate within the clade
	n.tip <- cache$n.tip # total number of tips in the tree
	len <- cache$len # branch length of each branch
	#init is a function that generates y=all possible scenarios of the tip initial states in a clade, target=idx of tip branches, t=tip branch length
	init <- function (group2, n.tip, len, nx, npad, states, states.sd) {
		target <- group2[group2<=n.tip] # find the tip branches
		# the initial state is y, which is a matrix with rows are tips, columns are c(E,DI,DR)
		if (length(target)==1) {
			y <- c(rep(0, nx),
                   rep(0, nx),
                   dnorm(x, states[target], states.sd[target]), rep(0, npad))
            list(target=target,t=len[target],y=y)
        } else if (length(target)>1) {
			y <- mapply(function(mean, sd) 
                          c(rep(0, nx),
                          dnorm(x, mean, sd), rep(0, npad),
                          rep(0, nx)),
                        states[target], states.sd[target], SIMPLIFY=FALSE)
            y[[1]][(2*nx+1):(3*nx)] <- y[[1]][(nx+1):(2*nx)]
            y[[1]][(nx+1):(2*nx)] <- 0
            list(target=target,t=len[target],y=y)
		} else {
			list(NULL)
		}
	}
	lapply(group,init,n.tip,len, nx, npad, states, states.sd) #apply init to each cell in the group list.
}

make.initial.conditions.protracted.mutation <- function(control) {
  tc <- control$tc
  r <- control$r
  nx.lo <- control$nx
  nx.hi <- nx.lo * r

  ## There is the chance that we could be slightly off on the depth
  ## by rounding error.  Because of this, I've done the testing
  ## against the *length* of the data, and then checked that the time
  ## is appropriate (to within eps of the correct value).  It is
  ## possible that two different branches with different numbers of
  ## nodes that finish right at the critical interval might have
  ## conflicting lengths.
  eps <- 1e-8
  function(init, pars, t, idx) {
    if ( length(init[[1]]) != length(init[[2]]) )
      stop("Data have incompatible length")

    if ( t < tc ) {
      ## if ( length(init[[1]]) / 2 == nx.hi ) { # t < tc
      ## if ( !((t - eps) < tc) )
      ##   stop("Wrong data size")
      nx <- nx.hi
      b <- pars[[1]]$b
    } else {
      ## if ( !((t + eps) > tc) )
      ##   stop("Wrong data size")
      nx <- nx.lo
      b <- pars[[2]]$b
    }
    ndat <- length(b)
    i <- seq_len(nx)
    j <- seq.int(nx+1, nx + ndat)
    z <- seq.int(nx+ndat+1, nx + 2*ndat)

    c(init[[1]][i],
      init[[1]][j] * init[[2]][j] * b, rep.int(0.0, nx - ndat), 
      (init[[1]][j] * init[[2]][z] + init[[2]][j] * init[[1]][z]) * b, rep.int(0.0, nx - ndat))
  }
}

## 5: root func
make.rootfunc.protracted.mutation <- function(cache) {
  root.idx <- cache$root
  nx <- cache$control$nx
  dx <- cache$control$dx
  
  function(res, pars, condition.surv, root, root.f, intermediates) {
    vals <- matrix(res$vals, nx, 3)[seq_len(pars$lo$ndat),]
    lq <- res$lq

    d.root <- vals[,3]

    root.p <- root.p.quasse(d.root, pars$lo, root, root.f)
    if ( condition.surv ) {
      b <- pars$lo$b
      e.root <- vals[,1]
      d.root <- d.root / sum(root.p * b * (1 - e.root)^2) * dx
    }

    loglik <- log(sum(root.p * d.root) * dx) + sum(lq)

    loglik
  }
}

######################################################################
## Extra core stuff:
#this is similar to make.all.branches.dtlik
make.all.branches.protracted.mutation <- function (cache,control) {
	branches <- make.branches.protracted.mutation(cache,control)
	initial.conditions <- make.initial.conditions.protracted.mutation(control)
	function (pars) {
		cache$y <- initial.tip.protracted.mutation(cache, cache$control, pars[[1]]$x)
		all.branches.matrix.protracted(pars,cache,initial.conditions,branches)
	}
}

make.branches.protracted.mutation <- function(cache,control) {
	if ( control$method == "fftC" ) {
    	branches <- make.branches.protracted.mutation.fftC(control)
  	} else if ( control$method == "fftR" ) {
    	branches <- make.branches.protracted.mutation.fftR(control)
	}
}

check.pars.protracted.mutation <- function(b.x, mu.x, lambda.x, drift, diffusion) {
  if ( any(!is.finite(c(b.x, mu.x, lambda.x, drift, diffusion))) )
    stop("Non-finite/NA parameters")
  if ( any(b.x < 0) || any(mu.x < 0) || any(lambda.x < 0) || diffusion <= 0 )
    stop("Illegal negative parameters")
  if ( !any(b.x > 0) )
    stop("No positive lambda; cannot compute likelihood")
}