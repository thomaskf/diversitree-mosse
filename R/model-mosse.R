#types: nucleotide types of tip species
#ntypes: number of nucleotide types
#states: mean substitution rates of tip species
#states.sd: standard deviation of substitution rates of tip species
#lambda: function of speciation rate against substitution rates
#mu: extinction rate
#Q: transition rate matrix among nucleotide types
#root.f: root nucleotype frequencies at equilibirum
#dt.max: define time interval, fix the decimal of branch length, so that we only need to calculate exp(xQdt) once. default is 0.01.
#xmax: the order of the upper limit of substitution rate. This is set to 10^-(order+3)*nx*r, so that each interval dx is fixed as 10^-(order+3), which make matrix exponential faster. The lower limit is 0. 
#nx: the number of bins for substituion rate, should be an integer power of 2. default is 1024.

make.mosse <- function (tree, types, ntypes=4, states, states.sd, lambda, mu, Q, control=NULL) {
	cache <- make.cache.mosse(tree, types, ntypes, states, states.sd, lambda, mu, Q, control)
	all.branches <- make.all.branches.mosse(cache,cache$control)
	rootfunc <- make.rootfunc.mosse(cache)
	f.pars <- make.pars.mosse(cache)
	ll <- function (pars, condition.surv=TRUE, root=ROOT.EQUI, root.f=NULL) {
		pars2 <- f.pars(pars)
		ans <- all.branches(pars2)
		rootfunc(ans, pars2, condition.surv, root, root.f, ntypes)
	}
	class(ll) <- c("mosse","dtlik","function")
	ll
}

## 2: info
make.info.mosse <- function(lambda, mu, phy) {
  ## Work around for .split:
  if ( !is.null(lambda) )
    argnames <- default.argnames.quasse(lambda, mu)
  else
    argnames <- NULL
  list(name="mosse",
       name.pretty="MoSSE",
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
       reference=NULL,
       ## These are special to QuaSSE:
       lambda=lambda,
       mu=mu)
}

## 3: make.cache
make.cache.mosse <- function(tree, types, ntypes, states, states.sd, lambda, mu, Q, control) {
  ## 1: tree
  #tree <- check.tree(tree)  

  ## 2: states & errors
  tmp <- check.states.quasse(tree, states, states.sd)
  states <- tmp$states
  states.sd <- tmp$states.sd

  ## 3: Control structure (lots of checking!)
  control <- check.control.mosse(control, tree, states)

  cache <- make.cache(tree)
  cache$edge.length <- round(cache$edge.length,digits=control$xmax+3)
  cache$depth <- round(cache$depth,digits=control$xmax+3)
  cache$types <- types
  cache$states  <- states
  cache$states.sd <- states.sd
  control$ntypes <- ntypes
  cache$control <- control

  ## 4: Speciation/extinction functions
  n.lambda <- check.f.quasse(lambda)
  n.mu     <- check.f.quasse(mu)
  n.args   <- n.lambda + n.mu + 2
  args <- list(lambda=seq_len(n.lambda),
               mu=seq_len(n.mu) + n.lambda,
               drift=n.lambda + n.mu + 1,
               diffusion=n.lambda + n.mu + 2)

  cache$lambda <- lambda
  cache$mu <- mu
  cache$args <- args
  cache$Q <- Q
  cache$info <- make.info.mosse(lambda, mu, tree)

  cache
}

## 4: initial.conditions
initial.tip.mosse <- function(cache, control, x) {
  nx <- control$nx * control$r
  npad <- nx - length(x)
  n <- control$ntypes
  
  message("n = ", n)
  message("nx = ", nx)
  message("npad = ", npad)
  message("x's length = ", length(x))
  message("x = ", capture.output(str(x)))
  message("cache$types = ", capture.output(cache$types))
  message("cache$states = ", capture.output(cache$states))
  message("cache$states.sd = ", capture.output(cache$states.sd))
  message("n = ", n)
  
  init <- function(type, mean, sd, n) {
    xlen <- length(x)
    if (xlen < 2) stop("x must have at least 2 points")
    out <- rep(0, nx * (n + 1))
    # extend x by one step so we have xlen+1 breakpoints -> xlen intervals
    dx_last <- x[xlen] - x[xlen - 1]
    x_ext <- c(x, x[xlen] + dx_last)
    # probability for each [x[i], x[i+1]) interval
    probs <- diff(pnorm(x_ext, mean, sd))  # length xlen
    out[nx*type+1:nx] <- c(probs, rep(0, npad))
    out
  }

    y <- mapply(init, cache$types, cache$states, cache$states.sd, n, SIMPLIFY=FALSE)
  dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
  # message("y$A's length: ", capture.output(length(y$A)))
  # message("y$B's length: ", capture.output(length(y$B)))
  # message("y$A[4097:4107]: ", capture.output(y$A[4097:4107]))
  # message("y$A[6264:6265]: ", capture.output(y$A[6264:6265]))
  # message("y$B[8197:8207]: ", capture.output(y$B[8197:8207]))
  # message("y$B[9860:9861]: ", capture.output(y$B[9860:9861]))
}

make.initial.conditions.mosse <- function(control) {
  tc <- control$tc
  r <- control$r
  nx.lo <- control$nx
  nx.hi <- nx.lo * r
  ntypes <- control$ntypes

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
      lambda <- pars[[1]]$lambda
    } else {
      ## if ( !((t + eps) > tc) )
      ##   stop("Wrong data size")
      nx <- nx.lo
      lambda <- pars[[2]]$lambda
    }
    
    ndat <- length(lambda)
    i <- seq_len(nx)
    j <- seq_len(ndat)
    out <- numeric(nx*(ntypes+1))
    out[i] <- init[[1]][i]
    for (z in 1:ntypes) {
    	out[nx*z+j] <- init[[1]][nx*z+j] * init[[2]][nx*z+j] * lambda
    }
    out
   }
}
   
## 5 rootfunc
## This function assumes that the root node is in the low-condition,
## which is enforced by the checking.
make.rootfunc.mosse <- function(cache) {
  root.idx <- cache$root
  nx <- cache$control$nx * cache$control$r
  dx <- cache$control$dx 
  # nx <- cache$control$nx
  # dx <- cache$control$dx * cache$control$r
  ntypes <- cache$control$ntypes
  Q <- cache$Q
  
  function(res, pars, condition.surv, root, root.f, ntypes) {

    # change to high (for debugging purpose)
    vals <- matrix(res$vals, nx, (ntypes+1))[seq_len(pars$hi$ndat),]
    # vals <- matrix(res$vals, nx, (ntypes+1))[seq_len(pars$lo$ndat),]
    
    message("vals[1,]", capture.output(vals[1,]));
    message("vals[2,]", capture.output(vals[2,]));
    message("vals[3,]", capture.output(vals[3,]));
    message("vals[4,]", capture.output(vals[4,]));
    message("vals[5,]", capture.output(vals[5,]));
    
    lq <- res$lq

    d.root <- vals[,-1]
    
    message("d.root[1,]", capture.output(d.root[1,]));
    message("d.root[2,]", capture.output(d.root[2,]));
    message("d.root[3,]", capture.output(d.root[3,]));
    message("d.root[4,]", capture.output(d.root[4,]));
    message("d.root[5,]", capture.output(d.root[5,]));
    
    # change to high (for debugging purpose)
	  # root.p <- root.p.mosse(d.root, pars$lo, root, root.f, ntypes)
    root.p <- root.p.mosse(d.root, pars$hi, root, root.f, ntypes)
    
    message("condition.surv:", capture.output(condition.surv));
    
    if ( condition.surv ) {
      # change to high (for debugging purpose)
      lambda <- pars$hi$lambda
      # lambda <- pars$lo$lambda
      e.root <- vals[,1]
      message("e.root", capture.output(e.root[1:6]));
      message("lambda's length", capture.output(length(lambda)))
      d.root <- mapply(function (i) d.root[,i] / (lambda * (1 - e.root)^2), 1:ntypes)
    }
    
    message("2 d.root[1,]", capture.output(d.root[1,]));
    message("2 d.root[2,]", capture.output(d.root[2,]));
    message("2 d.root[3,]", capture.output(d.root[3,]));
    message("2 d.root[4,]", capture.output(d.root[4,]));
    message("2 d.root[5,]", capture.output(d.root[5,]));
    
    message("dx", capture.output(dx));
    a <- log(sum(root.p * d.root) * dx)
    message("a", capture.output(a));
    message("lq", capture.output(lq));

        log(sum(root.p * d.root) * dx) + sum(lq)
  }
}

root.p.mosse <- function(d.root, pars, root, root.f, ntypes) {
  if ( !is.null(root.f) && root != ROOT.GIVEN )
    warning("Ignoring specified root state")
  
  x <- pars$x
  dx <- x[2] - x[1]
  
  message("pars$nx:", capture.output(pars$nx));
  message("ntypes:", capture.output(ntypes));
  message("dx:", capture.output(dx));
  message("d.root:", capture.output(d.root));

  if ( root == ROOT.FLAT ) {
    message("root == ROOT.FLAT");
    p <- 1 / ((pars$nx-1) * ntypes * dx)
  } else if ( root == ROOT.OBS )  {
    message("root == ROOT.OBS");
    p <- d.root / (sum(d.root) * dx)
  } else {
    root.i <- solve(t(cbind(c(1,1,1,1),pars$Q_orig[,-1])),c(1,0,0,0))
    message("root.i:", capture.output(root.i));
  	if ( root == ROOT.EQUI ) {
  	  message("root == ROOT.EQUI");
  	  p <- mapply(function (i) root.i[i] * d.root[,i] / (sum(d.root[,i]) * dx), 1:ntypes)
    } else if ( root == ROOT.GIVEN ){
      message("root == ROOT.GIVEN");
      p <- mapply(function (i) root.i[i] * root.f(x), 1:ntypes)
    }    	
  }
  
  message("p", capture.output(p));
  
  p
}

##6. make.pars
make.pars.mosse <- function(cache) {
  args <- cache$args
  dt.max <- cache$control$dt.max
  Q <- cache$Q
  ntypes <- cache$control$ntypes
  dx <- cache$control$dx
  dt.max <- cache$control$dt.max  
  r <- cache$control$r  
  n.hi <- cache$control$nx * r
  
  function(pars) {
    names(pars) <- NULL # Because of use of do.call, strip names

    drift <- pars[args$drift]
    diffusion <- pars[args$diffusion]

    ext <- mosse.extent(cache$control, drift, diffusion)

    ## Parameters, expanded onto the extent:
    pars <- expand.pars.quasse(cache$lambda, cache$mu, args, ext, pars)

    check.pars.quasse(pars$hi$lambda, pars$hi$mu, drift, diffusion)
	
	ndat <- pars$hi$ndat
	padding <- pars$hi$padding
	
	pars$hi$Q <- vector("list",n.hi-padding[2])
	pars$hi$Q[[1]] <- as.matrix(Matrix::expm(dx*Q*dt.max))
	for (i in 1:(n.hi-padding[2]-1)) {
		pars$hi$Q[[i+1]] <- pars$hi$Q[[i]]%*%pars$hi$Q[[1]]
	}
	pars$hi$Q <- pars$hi$Q[(padding[1]+1:ndat)]
	pars$lo$Q <- pars$hi$Q[seq(r,n.hi-padding[2],r)]
	
	if (cache$control$method == "fftC") {
		pars$hi$Q <- unlist(pars$hi$Q)
		pars$lo$Q <- unlist(pars$lo$Q)
	}
	pars$hi$Q_orig <- pars$lo$Q_orig <- Q
	
	pars
  }
}

mosse.extent <- function(control, drift, diffusion) {
  nx <- control$nx
  dx <- control$dx
  dt <- control$dt.max
  r  <- control$r
  w <- control$w

  mean <- drift * dt
  sd   <- sqrt(diffusion * dt)

    ## Another option here is to compute all the possible x values and
    ## then just drop the ones that are uninteresting?
  nkl <- ceiling(-(mean - w * sd)/dx/r) * c(r,1)
  nkr <- ceiling( (mean + w * sd)/dx/r) * c(r,1)
  ndat <- nx*c(r, 1) - (nkl + 1 + nkr)

  padding <- cbind(nkl, nkr)
  storage.mode(padding) <- "integer"

  ## Concatenate the x values, so that the lambda(x), mu(x)
  ## calculations work for both spaces simultaneously.
  x <- list(seq(dx, length.out=nx*r, by=dx)[(nkl[1]+1:ndat[1])],
            seq(dx*r, length.out=nx, by=dx*r)[(nkl[2]+1:ndat[2])])

  tr <- seq(r, length.out=ndat[2], by=r)

  list(x=x, padding=padding, ndat=ndat, tr=tr, nx=c(nx*r, nx))
}

##check control
check.control.mosse <- function(control, tree, states) {
  tree.length <- max(branching.times(tree))
  defaults <- list(tc=0.5,  #tree.length/10,
  				   dt.max=0.01,
                   nx=1024,
                   dx=10^(-(control$xmax+3)),
                   r=4,
                   w=5,
                   method="fftC",
                   flags=FFTW.MEASURE, # fftC only
                   atol=1e-6, # mol only
                   rtol=1e-6, # mol only
                   eps=1e-6,  # perhaps scale with dx?
                   verbose=FALSE)
  control <- if ( is.null(control) )
    defaults else modifyList(defaults, control)

  ## Eventually, this will contain "mol"
  method <- match.arg(control$method, c("fftC", "fftR"))

  # if ( control$tc <= 0 || control$tc >= tree.length )
  #  stop(sprintf("tc must lie in (0, %2.2f)", tree.length))
  if ( log2(control$nx) %% 1 != 0 )
    stop("nx must be a power of two")
  if ( log2(control$r) %% 1 != 0 )
    stop("r must be a power of two")

  ## These will be passed through to some C code, so type safety is
  ## important.
  ctrl.int <- c("nx", "flags", "verbose")
  ctrl.num <- c("tc", "dt.max", "r", "w", "atol", "rtol")
  control[ctrl.int] <- sapply(control[ctrl.int], as.integer)
  control[ctrl.num] <- sapply(control[ctrl.num], as.numeric)

  control
}

######################################################################
##core stuff:
make.all.branches.mosse <- function(cache, control) {
  branches <- make.branches.mosse(cache, control)
  initial.conditions <- make.initial.conditions.mosse(control)
  function(pars, preset=NULL) {
    cache$y <- initial.tip.mosse(cache, cache$control, pars[[1]]$x)
    all.branches.list(pars, cache, initial.conditions,
                      branches, preset)
  }
}

make.branches.mosse <- function(cache, control) {
  ## TODO: strictly, these should be backends...
  if ( control$method == "fftC" )
    branches <- make.branches.mosse.fftC(control)
  else if ( control$method == "fftR" )
    branches <- make.branches.mosse.fftR(control)
  else # already checked.
    stop("Unknown method", control$method)
}
