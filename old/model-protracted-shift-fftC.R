FFTW.ESTIMATE <- -1
FFTW.MEASURE <- 0
FFTW.PATIENT <- 1
FFTW.EXHAUSTIVE <- 2

## 8: branches
## TODO: I am not sure that this will pick up on the edge case of tc
## being *exactly* t0 or t0 + len (and what happens if rounding error
## puts this around a break).
## This might be something to check during the checking phase that tc
## does not come within 1e-8 of a node?
make.branches.protracted.shift.fftC <- function(control) {
  check.fftC()
  
  nx <- control$nx
  dx <- control$dx
  r <- control$r
  dt.max <- control$dt.max
  tc <- control$tc
  flags <- control$flags
  h <- control$h
  
  f.hi <- make.pde.protracted.shift.fftC(nx*r, dx/r, dt.max, 2*h+1, flags)
  f.lo <- make.pde.protracted.shift.fftC(nx,   dx,   dt.max, 2*h+1, flags)
  combine.branches.protracted(f.hi, f.lo, control)
}

make.pde.protracted.shift.fftC <- function(nx, dx, dt.max, nd, flags) {
  ptr <- .Call(r_make_protracted_shift_fft, as.integer(nx), as.numeric(dx),
               as.integer(nd), as.integer(flags))
  function(y, len, pars, t0, star, dt=dt.max) {
    nt <- as.integer(ceiling(len / dt.max))
    dt <- len / nt
    if ( !(length(y) %in% (nd * nx)) )
      stop("Wrong size y")
    if ( any(is.na(y)) )
      stop("Illegal NA values in input")
    if ( length(pars$lambda) != length(pars$mu) ||
         length(pars$lambda) > (nx-3) )
      stop("Incorrect length pars")
    if ( pars$diffusion <= 0 )
      stop("Invalid diffusion parameter")
    
    ans <- .Call(r_do_integrate_protracted_shift,
                 ptr, y, as.integer(star), pars$b, pars$mu, pars$lambda,
                 pars$s, pars$drift, pars$diffusion,
                 nt, dt, pars$padding)

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

