## MCMC for protracted, using slice sampler for pars and gibbs sampler for node state

mcmc.protracted <- function(tree, states, states.sd, b, lambda, mu, control=NULL, x.init, edge.args, nsteps, w, prior=NULL, edge.prior,
                         sampler=sampler.slice.protracted, fail.value=-Inf,
                         lower=-Inf, upper=Inf, print.every=1,
                         save.file, save.every=0, save.every.dt=NULL,
                         previous=NULL, previous.tol=1e-4,
                         keep.func=TRUE,
                         ...) {                 	
  if ( is.null(sampler) )
    sampler <- sampler.slice.protracted
  
  if ( save.every > 0 || !is.null(save.every.dt) ) {
    if ( missing(save.file) )
      stop("save.file must be given if save.every > 0")
    save.type <- tools::file_ext(save.file)
    save.type <- match.arg(tolower(save.type), c("csv", "rds"))
    save.fun <- switch(save.type,
                       rds=saveRDS,
                       csv=function(x, f) write.csv(x, f, row.names=FALSE))
    save.file.bak <- paste(save.file, ".bak", sep="")
  }

  lik <- make.protracted(tree, states, states.sd, b, lambda, mu, control=NULL)
  
  n.previous <- if ( is.null(previous) ) 0 else nrow(previous)
  npar <- length(argnames(lik))
  nedge <- length(edge.args)
  if ( !is.null(previous) ) {
    if ( !inherits(previous, "mcmcsamples") )
      stop("Currently only mcmcsamples objects can be continued")
    if ( n.previous >= nsteps ) {
      warning("Chain already complete")
      return(previous)
    }
    if ( !is.null(x.init) )
      stop("x.init must be NULL if continuing")
    hist.pars <- matrix(NA, ncol=npar, nrow=nsteps)
    hist.edge <- matrix(NA, ncol=nedge, nrow=nsteps)
    hist.prob <- rep(NA, nsteps)
    
    hist.pars[seq_len(n.previous),] <- previous[,par.args]
    hist.pars[seq_len(n.previous),] <- previous[,edge.args]
    hist.prob[seq_len(n.previous)]  <- previous$p

    x.init <- hist.pars[n.previous,]
    edge.init <- hist.pars[n.previous,]
    y.prev <- hist.prob[n.previous]
  } else {
    hist.pars <- matrix(NA, ncol=npar, nrow=nsteps)
    hist.edge <- matrix(NA, ncol=nedge, nrow=nsteps)
    hist.prob <- rep(NA, nsteps)
    
    edge.init <- tree$edge.spec[edge.args]
  }
  
  if ( is.null(names(x.init)) )
    try(colnames(hist.pars) <- names(x.init) <- argnames(lik),
        silent=TRUE)
  else
    colnames(hist.pars) <- names(x.init)
  
  if ( is.null(names(edge.init)) )
  	colnames(hist.edge) <- names(edge.init) <- edge.args
    
  if ( is.null(prior) && is.null(edge.prior) ) {
	  	  posterior <- protect(function(x) lik(x),
                         fail.value.default=fail.value)
  } else if ( is.null(prior) && !is.null(edge.prior) ) {
	 	  posterior <- protect(function(x, edge) lik(x) + sum(log(ifelse(edge==1, edge.prior, 1-edge.prior))),
                         fail.value.default=fail.value) 
  } else if ( !is.null(prior) && is.null(edge.prior) ) {
	 	  posterior <- protect(function(x, edge) lik(x) +prior(x),
                         fail.value.default=fail.value) 
  } else {
	 	  posterior <- protect(function(x, edge) lik(x) + prior(x) + sum(log(ifelse(edge==1, edge.prior, 1-edge.prior))),
                         fail.value.default=fail.value)
  }                       

  lower <- check.par.length(lower, npar)
  upper <- check.par.length(upper, npar)
  w     <- check.par.length(w,     npar)

  check.bounds(lower, upper, x.init)

  y.init <- posterior(x=x.init, edge=edge.init, fail.value=NULL)
  if ( !is.finite(y.init) ||
      (!is.null(fail.value) && y.init == fail.value) )
    stop("Starting point must have finite probability")
  if ( n.previous > 0 && abs(y.prev - y.init) > previous.tol ) {
    msg <- paste("Cannot continue the chain:\n",
                 sprintf("\texpected posterior = %2.7, got %2.7f",
                         previous$p[n.previous], y.init))
    stop(msg)
  }

  class.str <- c(sprintf("mcmcsamples.%s", get.info(lik)$name),
                 "mcmcsamples", "data.frame")

  clean.hist <- function(pars, p) {
    out <- data.frame(i=seq_along(p), pars, p)
    class(out) <- class.str
    out
  }

  we.should.print <- make.every.so.often(print.every)
  we.should.save <- make.every.so.often(save.every, save.every.dt)

  mcmc.loop <- function() {
    for ( i in seq(n.previous+1, nsteps, by=1) ) {
      tmp <- sampler(posterior, x.init, edge.init, y.init, w,
                     lower, upper, control)
      x.init <- hist.pars[i,] <- tmp[[1]]
      y.init <- tmp[[2]]
      for (j in seq_along(edge.init)) {
      	edge.init.new <- edge.init
      	edge.init.new[j] <- 1 - edge.init.new[j]
      	tree.new <- tree
      	tree.new$edge.spec[edge.args] <- edge.init.new
      	lik.new <- make.protracted(tree.new, states, states.sd, b, lambda, mu, control=NULL)
    
	  	if ( is.null(prior) && is.null(edge.prior) ) {
	  	  posterior.new <- protect(function(x, edge) lik.new(x),
                         fail.value.default=fail.value)
	 	} else if ( is.null(prior) && !is.null(edge.prior) ) {
	 	  posterior.new <- protect(function(x, edge) lik.new(x) + sum(log(ifelse(edge.init==1, edge.prior, 1-edge.prior))),
                         fail.value.default=fail.value) 
        } else if ( !is.null(prior) && is.null(edge.prior) ) {
	 	  posterior.new <- protect(function(x, edge) lik.new(x) + prior(x),
                         fail.value.default=fail.value) 
        } else {
	 	  posterior.new <- protect(function(x, edge) lik.new(x) + prior(x) + sum(log(ifelse(edge.init==1, edge.prior, 1-edge.prior))),
                         fail.value.default=fail.value)
        }
                         
    	 y.init.new <- posterior.new(x.init, edge.init.new, fail.value=NULL)
    	 if (tree.new$edge.spec[edge.args[j]]==1) {
    	 	if (!is.null(edge.prior)) {
    	 		p.j <- y.init.new * edge.prior[j]
    	 		p.j <- p.j / (p.j + y.init * (1 - edge.prior[j]))
    	 	} else {
    	 		p.j <- y.init.new
    	 		p.j <- p.j / (p.j + y.init)
    	 	}
    	 } else {
    	 	if (!is.null(edge.prior)) {
    	 		p.j <- y.init * edge.prior[j]
    	 		p.j <- p.j / (p.j + y.init.new * (1 - edge.prior[j]))
    	 	} else {
    	 		p.j <- y.init
    	 		p.j <- p.j / (p.j + y.init.new)
    	 	}
    	 }
    	 edge.init[j] <- as.numeric(runif(1,0,1) <= p.j)
    	 if (edge.init[j]!=tree$edge.spec[edge.args[j]]) {
    	 	tree <- tree.new
    	 	lik <- lik.new
    	 	y.init <- y.init.new
    	 	
  if ( is.null(prior) && is.null(edge.prior) ) {
	  	  posterior <- protect(function(x) lik(x),
                         fail.value.default=fail.value)
  } else if ( is.null(prior) && !is.null(edge.prior) ) {
	 	  posterior <- protect(function(x, edge) lik(x) + sum(log(ifelse(edge==1, edge.prior, 1-edge.prior))),
                         fail.value.default=fail.value) 
  } else if ( !is.null(prior) && is.null(edge.prior) ) {
	 	  posterior <- protect(function(x, edge) lik(x) +prior(x),
                         fail.value.default=fail.value) 
  } else {
	 	  posterior <- protect(function(x, edge) lik(x) + prior(x) + sum(log(ifelse(edge==1, edge.prior, 1-edge.prior))),
                         fail.value.default=fail.value)
  }  
        
    	 }
      }
      hist.edge[i,] <- edge.init
	  hist.prob[i]  <- y.init
	  
      if ( we.should.print() )
          print(cat(i,
                    paste(sprintf("%2.4f", x.init), collapse=", "), paste(sprintf("%2.4f", edge.init), collapse=", "),
                    y.init))
      if ( we.should.save() ) {
        j <- seq_len(i)
        ## Back up the old version to avoid IO errors if the system
        ## fails while saving.
        if ( file.exists(save.file) )
          ok <- file.rename(save.file, save.file.bak)
        ok <- try(save.fun(clean.hist(hist.pars[j,,drop=FALSE], hist.prob[j]),
                           save.file))
        if ( inherits(ok, "try-error") )
          warning("Error while writing progress file (continuing)",
                  immediate.=TRUE)
      }
    }
    clean.hist(hist.pars, hist.edge, hist.prob)
  }

  mcmc.recover <- function(...) {
    j <- !is.na(hist.prob)
    if ( !any(j) )
      stop("MCMC was stopped before any samples taken")
    hist <- clean.hist(hist.pars[j,], hist.edge[j,], hist.prob[j])
    warning("MCMC was stopped prematurely: ", nrow(hist), "/", nsteps,
            " steps completed.  Truncated chain is being returned.",
            immediate.=TRUE)
    hist
  }

  samples <- tryCatch(mcmc.loop(), interrupt=mcmc.recover)

  if ( save.every > 0 || !is.null(save.every.dt) )
    if ( nrow(samples) == nsteps && file.exists(save.file.bak) )
      file.remove(save.file.bak)

  if (keep.func) {
    attr(samples, "func")  <- set.defaults(lik, defaults=list(...))
    attr(samples, "prior") <- prior
  }

  samples
}

sampler.slice.protracted <- function(lik, x.init, edge.init, y.init, w, lower, upper, control) {
  for ( i in seq_along(x.init) ) {
    xy <- slice.1d(make.unipar.protracted(lik, x.init, edge.init, i),
                   x.init[i], y.init, w[i], lower[i], upper[i])
    x.init[i] <- xy[1]
    y.init    <- xy[2]
  }

  list(x.init, y.init)
}

make.unipar.protracted <- function(f, x, edge, i) {
  function(z) {
    x[i] <- z
    f(x,edge)
  }
}

