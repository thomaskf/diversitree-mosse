starting.point.protracted.mutation <- function(tree, q.div=5, states, states.sd=NULL) {
  pars.bd <- starting.point.bd(tree)
  if  ( pars.bd[1] > pars.bd[2] )
    p <- c(pars.bd[1],(pars.bd[1] - pars.bd[2]) / q.div,pars.bd[1])
  else
    p <- c(pars.bd[1],pars.bd[1]/ q.div,pars.bd[1])
  lik.bm <- make.bm(tree, states, states.sd,
                    control=list(method="pruning", backend="C"))
  c(p, diffusion=as.numeric(stats::coef(find.mle(lik.bm, .1))))
}

## I use a number of elements of pars.
## pars[[i]]${b,mu,lambda,drift,diffusion,padding}
## pars$tr
expand.pars.protracted.mutation <- function(b, mu, lambda, args, ext, pars) {
  pars.use <- vector("list", 2)
  for ( i in c(1,2) ) {
    x <- list()
    pars.use[[i]] <-
      list(x=ext$x[[i]], # May screw other things up (was $x[i])
           b=do.call(b, c(ext$x[i], pars[args$b])),
           mu=do.call(mu, c(ext$x[i], pars[args$mu])),
           lambda=do.call(lambda, c(ext$x[i], pars[args$lambda])),
           drift=pars[args$drift],
           diffusion=pars[args$diffusion],
           padding=ext$padding[i,],
           ndat=ext$ndat[i],
           nx=ext$nx[i])
  }
  names(pars.use) <- c("hi", "lo")
  pars.use$tr <- ext$tr
  pars.use
}

make.pars.protracted.mutation <- function(cache) {
  args <- cache$args

  function(pars) {
    names(pars) <- NULL # Because of use of do.call, strip names

    drift <- pars[args$drift]
    diffusion <- pars[args$diffusion]

    ext <- quasse.extent(cache$control, drift, diffusion)
    ## This would confirm the translation:
    ##   all.equal(ext$x[[1]][ext$tr], ext$x[[2]])

    ## Parameters, expanded onto the extent:
    pars <- expand.pars.protracted.mutation(cache$b, cache$mu, cache$lambda, args, ext, pars)

    check.pars.protracted.mutation(pars$hi$b, pars$hi$mu, pars$hi$lambda, drift, diffusion)

    pars
  }
}

combine.branches.protracted <- function(f.hi, f.lo, control) {
  nx <- control$nx
  dx <- control$dx
  tc <- control$tc
  r <- control$r
  eps <- log(control$eps)
  dt.max <- control$dt.max

  ## This is hacky version of the log compensation.  It also reduces
  ## the stepsize when bisecting a branch.  It doesn't seem to change
  ## much on 
  careful <- function(f, y, len, pars, t0, star, dt.max) {
    ans <- f(y, len, pars, t0, star)
    if ( ans[[1]] > eps ) { # OK
      ans
    } else {
      if ( control$method == "fftC" ||
           control$method == "fftR" )
        dt.max <- dt.max / 2 # Possibly needed
      len2 <- len/2
      ans1 <- Recall(f, y,         len2, pars, t0,        dt.max)
      ans2 <- Recall(f, ans1[[2]], len2, pars, t0 + len2, dt.max)
      ans2[[1]][[1]] <- ans1[[1]][[1]] + ans2[[1]][[1]]
      ans2
    }
  }

  ## Start by normalising the input so that eps up there make
  ## sense...
  function(y, len, pars, t0, star) {
    if ( t0 < tc ) {
      dx0 <- dx / r
      nx0 <- nx * r
    } else {
      dx0 <- dx
      nx0 <- nx
    }

    ## Here, we also squash all negative numbers.
    if ( any(y < -1e-8) )
      stop("Actual negative D value detected -- calculation failure")
    y[y < 0] <- 0
    y <- matrix(y, nx0, 3)
    q0 <- sum(y[,c(2,3)]) * dx0
    if ( q0 <= 0 )
      stop("No positive D values")
    y[,c(2,3)] <- y[,c(2,3)] / q0
    lq0 <- log(q0)
    
    if ( t0 >= tc ) {
      ans <- careful(f.lo, y, len, pars$lo, t0, star, dt.max)
    } else if ( t0 + len < tc ) {
      ans <- careful(f.hi, y, len, pars$hi, t0, star, dt.max)
    } else {
      len.hi <- tc - t0
      ans.hi <- careful(f.hi, y, len.hi, pars$hi, t0, star, dt.max)

      y.lo <- ans.hi[[2]][pars$tr,]
      lq0 <- lq0 + ans.hi[[1]]
      if ( nrow(y.lo) < nx )
        y.lo <- rbind(y.lo, matrix(0, nx - length(pars$tr), 3))

      ## Fininshing up with the low resolution branch...
      ans <- careful(f.lo, y.lo, len - len.hi, pars$lo, tc, dt.max)
    }

    c(ans[[1]] + lq0, ans[[2]])
  }
}
