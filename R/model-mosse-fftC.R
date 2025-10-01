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
    
    message("length(y):", length(y));
    message("y[1:2,1]:", capture.output(y[1:2,1]));
    message("y[1:2,2]:", capture.output(y[1:2,2]));
    message("y[1:2,3]:", capture.output(y[1:2,3]));
    message("y[1:2,4]:", capture.output(y[1:2,4]));
    message("y[1:2,5]:", capture.output(y[1:2,5]));
    message("length(pars$lambda):", length(pars$lambda));
    message("pars$lambda:", capture.output(pars$lambda));
    message("length(pars$mu):", length(pars$mu));
    message("pars$mu:", capture.output(pars$mu));
    message("pars$drift:", capture.output(pars$drift));
    message("pars$diffusion:", capture.output(pars$diffusion));
    message("length(pars$Q):", length(pars$Q));
    message("pars$Q:", capture.output(pars$Q));
    message("nt:", capture.output(nt));
    message("dt.max:", capture.output(dt.max));
    message("pars$padding:", capture.output(pars$padding));
    
    ans <- .Call(r_do_integrate_mosse,
                 ptr, y, pars$lambda, pars$mu,
                 pars$drift, pars$diffusion, pars$Q, 
                 as.integer(nt), dt.max, as.integer(pars$padding))
    
    message("length(ans):", length(ans));
    message("ans[1:3]", capture.output(ans[1:3]));
    message("ans[4097:4099", capture.output(ans[4097:4099]));
    message("ans[8193:8195]", capture.output(ans[8193:8195]));
    message("ans[12289:12291", capture.output(ans[12289:12291]));
    message("ans[16385:16387]", capture.output(ans[16385:16387]));
    
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
