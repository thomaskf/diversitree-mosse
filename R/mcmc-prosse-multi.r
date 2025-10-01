## MCMC for prosse.multi, using slice sampler for pars and gibbs sampler for node state

mcmc.prosse.multi <- function(lik, tree, species.name, unknown.tip, traits, states, states.sd, lambda, control=NULL, x.init, nsteps, w, prior=NULL,
                         sampler=sampler.slice, fail.value=-Inf,
                         lower=-Inf, upper=Inf, print.every=1,
                         save.file, save.every=0, save.every.dt=NULL,
                         previous=NULL, previous.tol=1e-4,
                         keep.func=TRUE,
                         ...) {                 	
  if ( is.null(sampler) )
    sampler <- sampler.slice
  
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

  n.previous <- if ( is.null(previous) ) 0 else nrow(previous)
  
  par.args <- c(argnames(lik), unknown.tip)
  n.par <- length(argnames(lik))
  npar <- n.par + length(unknown.tip)
  
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
    hist.prob <- rep(NA, nsteps)
    
    hist.pars[seq_len(n.previous),] <- previous[,par.args]
    hist.prob[seq_len(n.previous)]  <- previous$p

    x.init <- hist.pars[n.previous,]
    y.prev <- hist.prob[n.previous]
  } else {
    hist.pars <- matrix(NA, ncol=npar, nrow=nsteps)
    hist.prob <- rep(NA, nsteps)
  }
  
  if ( is.null(names(x.init)) ) {
    try(colnames(hist.pars) <- names(x.init) <- par.args,
        silent=TRUE)
  } else {
    colnames(hist.pars) <- names(x.init)
  }
      
  if ( is.null(prior) ) {
      posterior <- protect(function(x) lik(x),
                         fail.value.default=fail.value)
  } else {
      posterior <- protect(function(x) lik(x) + prior(x),
                         fail.value.default=fail.value)
  }
                         
  lower <- check.par.length(lower, n.par)
  upper <- check.par.length(upper, n.par)
  w     <- check.par.length(w,     n.par)

  check.bounds(lower, upper, x.init[1:n.par])

  y.init <- posterior(x=x.init[1:n.par], fail.value=NULL)
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
      tree$species[unknown.tip] <- x.init[unknown.tip]
      for (j in unknown.tip) {
          species <- unique(c(tree$species,species.name[j]))
          nspecies <- length(species)
          tmp <- numeric(nspecies)
          for (z in 1:nspecies) {
              tree$species[j] <- species[z]
              if(check.paraphyletic(tree)) {
                  tmp[z] <- fail.value
              } else {
                  lik <- make.prosse.multi(tree, traits, states, states.sd, lambda, control=NULL)
                  lik <- constrain(lik,drift~0)
                  tmp[z] <- posterior(x=x.init[1:n.par], fail.value=NULL)
              }
          }
          x.init[j] <- hist.pars[i,j] <- sample(x=species,size=1,prob=exp(tmp))
      }
      tree$species[unknown.tip] <- x.init[unknown.tip]
      lik <- make.prosse.multi(tree, traits, states, states.sd, lambda, control=NULL)
      lik <- constrain(lik,drift~0)
      tmp <- sampler(posterior, x.init[1:n.par], y.init, w,
                     lower, upper, control)
      x.init[1:n.par] <- hist.pars[i,1:n.par] <- tmp[[1]]
      y.init <- hist.prob[i]  <- tmp[[2]]
	
      if ( we.should.print() )
          print(cat(i, paste(sprintf("%2.4f", x.init), collapse=", "), y.init))
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
    clean.hist(hist.pars, hist.prob)
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

getMRCA <- function (phy, tip) {
    if (length(tip) < 2) {return(NULL)}
    Ntip <- length(phy$tip.label)
    rootnd <- Ntip + 1L
    pars <- integer(phy$Nnode)
    tnd <- tip
    done_v <- logical(Ntip + phy$Nnode)
    pvec <- integer(Ntip + phy$Nnode)
    pvec[phy$edge[, 2]] <- phy$edge[, 1]
    nd <- tnd[1]
    for (k in 1:phy$Nnode) {
        nd <- pvec[nd]
        pars[k] <- nd
        if (nd == rootnd)
            break
    }
    pars <- pars[1:k]
    mrcind <- integer(max(pars))
    mrcind[pars] <- 1:k
    mrcand <- pars[1]
    for (i in 2:length(tnd)) {
        cnd <- tnd[i]
        done <- done_v[cnd]
        while (!done) {
            done_v[cnd] <- TRUE
            cpar <- pvec[cnd]
            done <- done_v[cpar]
            if (cpar %in% pars) {
                if (cpar == rootnd)
                  return(rootnd)
                if (mrcind[cpar] > mrcind[mrcand])
                  mrcand <- cpar
                done_v[cpar] <- TRUE
                done <- TRUE
            }
            cnd <- cpar
        }
    }
    mrcand
}
keep.tip.fixed <- function (tree, tip) {
    des <- unique(tree$edge[sapply(tip,function (i) which(tree$edge[,2]==i)),1])
    node <- des
    des <- des[des!=(length(tree$tip.label)+1)]
    while(length(des)>1) {
        des <- unique(tree$edge[sapply(des,function (i) which(tree$edge[,2]==i)),1])
        node <- c(node,des)
        des <- des[des!=(length(tree$tip.label)+1)]
    }
    mrca <- getMRCA(tree,tip)
    node[node>=mrca]
}
check.paraphyletic <- function (tree) {
    species <- unique(tree$species)
    nspecies <- length(species)
    clades <- sapply(species,function (i) keep.tip.fixed(tree, which(tree$species==i)))
    tmp <- sapply(1:nspecies, function (i) sapply(clades, function (j) sum(is.element(clades[[i]],j))>0))
    tmp <- sapply(1:nspecies,function(i) sapply(c(1:nspecies)[-i], function (j) tmp[i,j]&&tmp[j,i]))
    sum(tmp)>0
}
