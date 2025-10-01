make.branches.protracted.shift.fftR <- function(control) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc
  h <- control$h

  f.hi <- make.pde.protracted.shift.fftR(nx*r, dx/r, dt.max, 2*h+1)
  f.lo <- make.pde.protracted.shift.fftR(nx,   dx,   dt.max, 2*h+1)
  combine.branches.protracted(f.hi, f.lo, control)
}

make.pde.protracted.shift.fftR <- function(nx, dx, dt.max, nd) {
  function(y, len, pars, t0, star) {
    padding <- pars$padding
    ndat <- length(pars$lambda)
    
    nt <- as.integer(ceiling(len / dt.max))
    dt <- len / nt
    if ( !(length(y) %in% (nd * nx)) )
      stop("Wrong size y")
    if ( length(pars$lambda) != length(pars$mu) ||
         length(pars$lambda) > (nx-3) )
      stop("Incorrect length pars")
    if ( pars$diffusion <= 0 )
      stop("Invalid diffusion parameter")

    if ( !is.matrix(y) )
      y <- matrix(y, nx, nd)

    ans <- protracted.integrate.shift.fftR(y, star, pars$b, pars$mu, pars$lambda, pars$s, pars$drift,
                                 pars$diffusion, nt, dt, nx, ndat, dx, nd,
                                 padding[1], padding[2])
    ## Do the log compensation here, to make the careful calcuations
    ## easier later.
    q <- sum(ans[,-1]) * dx
    ans[,-1] <- ans[,-1] / q
    list(log(q), ans)
  }
}

protracted.integrate.shift.fftR <- function(vars, star, b, mu, lambda, s, drift, diffusion,
                                  nstep, dt, nx, ndat, dx, nd, nkl, nkr) {
  kern <- fftR.make.kern(-dt * drift, sqrt(dt * diffusion),
                         nx, dx, nkl, nkr)
  fy <- fft(kern)  
  for ( i in seq_len(nstep) ) {
    vars <- fftR.propagate.t.protracted.shift(vars, star, b, mu, lambda, s, dt, ndat, nd)
    vars <- fftR.propagate.x(vars, nx, fy, nkl, nkr)
  }
  vars
}

fftR.propagate.t.protracted.shift <- function(vars, star, b, mu, lambda, s, dt, ndat, nd) {
  i <- seq_len(ndat)
  r <- b - mu
  z <- exp(dt * r)
  e0 <- vars[i,1]
  h <- (nd-1)/2
  vars[i,1] <- (mu + z*(e0 - 1)*mu - b*e0) /
    (mu + z*(e0 - 1)*b - b*e0)
  A <- (z * r * r)/(z * b - mu + (1-z)*b*e0)^2
  if (star==1) {
  	vars[i,-1] <- A * vars[i,-1] *exp(-lambda * dt)
  } else if (star==2) {
  	dr0 <- vars[i,(h+1):(2*h+1)]
  	sum.dr0 <- sum(dr0)/h
  	vars[i,(h+1):(2*h+1)] <- ((dr0-sum.dr0) * exp(-s*lambda*h * dt) + sum.dr0) * A
  	vars[i,2:(h+1)] <- vars[i,(h+1):(2*h+1)] + (vars[i,2:(h+1)] - dr0) * exp(-(s*lambda*(h-1)+lambda) * dt) * A
  }
  vars
}