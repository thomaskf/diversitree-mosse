starting.point.protracted.shift <- function(tree, q.div=5, states, states.sd=NULL) {
  pars.bd <- starting.point.bd(tree)
  if  ( pars.bd[1] > pars.bd[2] )
    p <- c(pars.bd[1],(pars.bd[1] - pars.bd[2]) / q.div,pars.bd[1],pars.bd[1])
  else
    p <- c(pars.bd[1],pars.bd[1]/ q.div,pars.bd[1],pars.bd[1])
  lik.bm <- make.bm(tree, states, states.sd,
                    control=list(method="pruning", backend="C"))
  c(p, diffusion=as.numeric(stats::coef(find.mle(lik.bm, .1))))
}

## I use a number of elements of pars.
## pars[[i]]${b,mu,lambda,drift,diffusion,padding}
## pars$tr
expand.pars.protracted.shift <- function(b, mu, lambda, args, ext, pars) {
  pars.use <- vector("list", 2)
  for ( i in c(1,2) ) {
    x <- list()
    pars.use[[i]] <-
      list(x=ext$x[[i]], # May screw other things up (was $x[i])
           b=do.call(b, c(ext$x[i], pars[args$b])),
           mu=do.call(mu, c(ext$x[i], pars[args$mu])),
           lambda=do.call(lambda, c(ext$x[i], pars[args$lambda])),
           s=pars[args$s],
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

make.pars.protracted.shift <- function(cache) {
  args <- cache$args

  function(pars) {
    names(pars) <- NULL # Because of use of do.call, strip names
	
	s <- pars[args$s]
    drift <- pars[args$drift]
    diffusion <- pars[args$diffusion]

    ext <- quasse.extent(cache$control, drift, diffusion)
    ## This would confirm the translation:
    ##   all.equal(ext$x[[1]][ext$tr], ext$x[[2]])

    ## Parameters, expanded onto the extent:
    pars <- expand.pars.protracted.shift(cache$b, cache$mu, cache$lambda, args, ext, pars)

    check.pars.protracted.shift(pars$hi$b, pars$hi$mu, pars$hi$lambda, s, drift, diffusion)

    pars
  }
}
