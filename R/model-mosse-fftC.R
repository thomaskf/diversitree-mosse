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
