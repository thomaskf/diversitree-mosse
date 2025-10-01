make.branches.protracted.mutation.fftR <- function(control) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc

  f.hi <- make.pde.protracted.mutation.fftR(nx*r, dx/r, dt.max, 3L)
  f.lo <- make.pde.protracted.mutation.fftR(nx,   dx,   dt.max, 3L)
  combine.branches.protracted(f.hi, f.lo, control)
}

make.pde.protracted.mutation.fftR <- function(nx, dx, dt.max, nd) {
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

    ans <- protracted.integrate.mutation.fftR(y, star, pars$b, pars$mu, pars$lambda, pars$drift,
                                 pars$diffusion, nt, dt, nx, ndat, dx,
                                 padding[1], padding[2])
    ## Do the log compensation here, to make the careful calcuations
    ## easier later.
    q <- sum(ans[,-1]) * dx
    ans[,-1] <- ans[,-1] / q
    list(log(q), ans)
  }
}

protracted.integrate.mutation.fftR <- function(vars, star, b, mu, lambda, drift, diffusion,
                                  nstep, dt, nx, ndat, dx, nkl, nkr) {
  kern <- fftR.make.kern(-dt * drift, sqrt(dt * diffusion),
                         nx, dx, nkl, nkr)
  fy <- fft(kern)  
  for ( i in seq_len(nstep) ) {
    vars <- fftR.propagate.t.protracted.mutation(vars, star, b, mu, lambda, dt, ndat)
    vars <- fftR.propagate.x(vars, nx, fy, nkl, nkr)
  }
  vars
}

fftR.propagate.t.protracted.mutation <- function(vars, star, b, mu, lambda, dt, ndat) {
  i <- seq_len(ndat)
  r <- b - mu
  z <- exp(dt * r)
  e0 <- vars[i,1]
  di0 <- vars[i,2]
  dr0 <- vars[i,3]
  vars[i,1] <- (mu + z*(e0 - 1)*mu - b*e0) /
    (mu + z*(e0 - 1)*b - b*e0)
  A <- (z * r * r)/(z * b - mu + (1-z)*b*e0)^2
  if (star==1) {
  	vars[i,2] <- A * di0 * exp(-lambda * dt)
  	vars[i,3] <- A * dr0 * exp(-lambda * dt)
  } else if (star==2) {
  	vars[i,3] <- A * dr0
  	vars[i,2] <- vars[i,3] + (di0-dr0) * A / exp(dt * (b + mu))
  }
  vars
}
