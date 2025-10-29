FFTW.ESTIMATE <- -1
FFTW.MEASURE <- 0
FFTW.PATIENT <- 1
FFTW.EXHAUSTIVE <- 2

make.branches.mosse.fftC <- function(control) {
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc
  ntypes <- control$ntypes
  flags <- control$flags

  f.hi <- make.pde.mosse.fftC(nx*r, dx, dt.max, ntypes+1, flags)
  f.lo <- make.pde.mosse.fftC(nx,   dx * r,   dt.max, ntypes+1, flags)
  combine.branches.mosse(f.hi, f.lo, control)
}

make.pde.mosse.fftC <- function(nx, dx, dt.max, nd, flags) {
  
  message(".")
  message("nx = ", capture.output(nx))
  message("dx = ", capture.output(dx))
  message("nd = ", capture.output(nd))
  message("flags = ", capture.output(flags))
  
  ptr <- .Call(r_make_mosse_fft, as.integer(nx), as.numeric(dx),
               as.integer(nd), as.integer(flags))
  function(y, len, pars, t0, dt.max) {
    nt <- as.integer(ceiling(len / dt.max))
    if ( !(length(y) %in% (nd * nx)) )
      stop("Wrong size y")
    if ( any(is.na(y)) )
      stop("Illegal NA values in input")
    if ( length(pars$lambda) != length(pars$mu) ||
         length(pars$lambda) > (nx-3) )
      stop("Incorrect length pars")
    if ( pars$diffusion <= 0 )
      stop("Invalid diffusion parameter")
    
    message(".")
    message("length(y) = ", capture.output(length(y)))
    x <- length(y)/5
    message("y[1:5] = ", capture.output(y[1:5]))
    s <- x+1; t <- x+5;
    message("y[", s, ":", t , "] = ", capture.output(y[s:t]))
    s <- 2*x+1; t <- 2*x+5;
    message("y[", s, ":", t , "] = ", capture.output(y[s:t]))
    s <- 3*x+1; t <- 3*x+5;
    message("y[", s, ":", t , "] = ", capture.output(y[s:t]))
    s <- 4*x+1; t <- 4*x+5;
    message("y[", s, ":", t , "] = ", capture.output(y[s:t]))
    message("length(pars$lambda) = ", capture.output(length(pars$lambda)))
    message("pars$lambda[1:5] = ", capture.output(pars$lambda[1:5]))
    message("length(pars$mu) = ", capture.output(length(pars$mu)))
    message("pars$mu[1:5] = ", capture.output(pars$mu[1:5]))
    message("pars$drift = ", capture.output(pars$drift))
    message("pars$diffusion = ", capture.output(pars$diffusion))
    message("length(pars$Q) = ", capture.output(length(pars$Q)))
    message("nt = ", capture.output(nt))
    message("dt.max = ", capture.output(dt.max))
    message("pars$padding = ", capture.output(pars$padding))
    
    ans <- .Call(r_do_integrate_mosse,
                 ptr, y, pars$lambda, pars$mu,
                 pars$drift, pars$diffusion, pars$Q, 
                 as.integer(nt), dt.max, as.integer(pars$padding))

    ## Do the log compensation here, to make the careful calcuations
    ## easier later.
    if ( ncol(ans) > 1 ) {
      q <- sum(ans[,-1]) * dx
      ans[,-1] <- ans[,-1] / q
      list(log(q), ans)
    } else {
      list(0, ans)
    }
  }
}

check.fftC <- function(error=TRUE) {
  ok <- is.loaded("r_make_mosse_fft", "diversitree")
  if ( error && !ok )
    stop("diversitree was built without FFTW support.  ",
         "See Details in ?make.quasse")
  ok
}
