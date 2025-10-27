make.branches.mosse.fftR <- function(control) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc
  ntypes <- control$ntypes

  f.hi <- make.pde.mosse.fftR(nx*r, dx, dt.max, ntypes+1L)
  f.lo <- make.pde.mosse.fftR(nx,   dx * r,   dt.max, ntypes+1L)
  combine.branches.mosse(f.hi, f.lo, control)
}

make.pde.mosse.fftR <- function(nx, dx, dt.max, nd) {
  function(y, len, pars, t0, dt.max) {
    padding <- pars$padding
    ndat <- length(pars$lambda)
    
    nt <- as.integer(ceiling(len / dt.max))
    if ( !(length(y) %in% (nd * nx)) )
      stop("Wrong size y")
    if ( length(pars$lambda) != length(pars$mu) ||
         length(pars$lambda) > (nx-3) )
      stop("Incorrect length pars")
    if ( pars$diffusion <= 0 )
      stop("Invalid diffusion parameter")

    if ( !is.matrix(y) )
      y <- matrix(y, nx, nd)

    ans <- mosse.integrate.fftR(y, pars$lambda, pars$mu, pars$drift,
                                 pars$diffusion, pars$Q, nt, dt.max, nx, ndat, dx,
                                 padding[1], padding[2])
    ## Do the log compensation here, to make the careful calcuations
    ## easier later.
    q <- sum(ans[,-1]) * dx
    ans[,-1] <- ans[,-1] / q
    list(log(q), ans)
  }
}

## Note that the sign of drift here needs inverting.  (see
## src/quasse-eqs-fftC.c:qf_setup_kern() for the equivalent flip).
mosse.integrate.fftR <- function(vars, lambda, mu, drift, diffusion, Q,
                                  nstep, dt, nx, ndat, dx, nkl, nkr) {
  kern <- fftR.make.kern(-dt * drift, sqrt(dt * diffusion),
                         nx, dx, nkl, nkr)
  fy <- fft(kern)  
  for ( i in seq_len(nstep) ) {
    vars <- fftR.mosse.propagate.t(vars, lambda, mu, Q, dt, ndat)
    vars <- fftR.propagate.x(vars, nx, fy, nkl, nkr)
  }
  vars
}

fftR.mosse.propagate.t <- function(vars, lambda, mu, Q, dt, ndat) {
  i <- seq_len(ndat)
  r <- lambda - mu
  z <- exp(dt * r)
  e0 <- vars[i,1]
  d0 <- vars[i,-1]
  vars[i,1] <- (mu + z*(e0 - 1)*mu - lambda*e0) /
    (mu + z*(e0 - 1)*lambda - lambda*e0)
  dd <- (z * r * r)/(z * lambda - mu + (1-z)*lambda*e0)^2
  d0 <- dd * d0
  vars[i,-1] <- t(mapply(function (z) d0[z,] %*% Q[[z]], i))
  vars
}

combine.branches.mosse <- function(f.hi, f.lo, control) {
  nx <- control$nx
  dx <- control$dx
  tc <- control$tc
  r <- control$r
  eps <- log(control$eps)
  dt.max <- control$dt.max
  nd <- control$ntypes+1

  ## Start by normalising the input so that eps up there make
  ## sense...
  function(y, len, pars, t0, idx) {
    if ( t0 < tc ) {
      dx0 <- dx
      nx0 <- nx * r
    } else {
      dx0 <- dx * r
      nx0 <- nx
    }

    ## Here, we also squash all negative numbers.
    if ( any(y < -1e-8) )
      stop("Actual negative D value detected -- calculation failure")
    y[y < 0] <- 0
    y <- matrix(y, nx0, length(y)/nx0)
    q0 <- sum(y[,-1]) * dx0
    if ( q0 <= 0 )
      stop("No positive D values")
    y[,-1] <- y[,-1] / q0
    lq0 <- log(q0)
    
    if ( t0 >= tc ) {
      ans <- f.lo(y, len, pars$lo, t0, dt.max)
    } else if ( t0 + len < tc ) {
      ans <- f.hi(y, len, pars$hi, t0, dt.max)
    } else {
      len.hi <- tc - t0
      ans.hi <- f.hi(y, len.hi, pars$hi, t0, dt.max)
      message("len.hi = ", capture.output(len.hi))
      message("length(pars$tr) = ", capture.output(length(pars$tr)))
      message("pars$tr[1:5] = ", capture.output(pars$tr[1:5]))
      
      y.lo <- ans.hi[[2]][pars$tr,]
      lq0 <- lq0 + ans.hi[[1]]
      if ( nrow(y.lo) < nx )
        y.lo <- rbind(y.lo, matrix(0, nx - length(pars$tr),nd))
        
      ## Fininshing up with the low resolution branch...
      message("len - len.hi = ", capture.output(len - len.hi))
      ans <- f.lo(y.lo, len - len.hi, pars$lo, tc, dt.max)
    }

    c(ans[[1]] + lq0, ans[[2]])
  }
}
