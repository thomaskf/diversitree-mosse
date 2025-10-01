#include "config.h"

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> /* why? */
#include <Rmath.h> /* for dnorm() */

#ifdef HAVE_FFTW3_H
#include <complex.h>
#include <fftw3.h>
#include "rfftw.h"
#include "quasse-eqs-fftC.h"
#include "mosse-eqs-fftC.h"

static void mosse_fft_finalize(SEXP extPtr);

SEXP r_make_mosse_fft(SEXP r_nx, SEXP r_dx, SEXP r_nd, SEXP r_flags) {
  mosse_fft *obj;
  SEXP extPtr;
  int nx = INTEGER(r_nx)[0];
  double dx = REAL(r_dx)[0];
  int n_fft = LENGTH(r_nd);
  int i;
  int flags = FFTW_MEASURE;
  int *nd = (int*)calloc(n_fft, sizeof(int));
  for ( i = 0; i < n_fft; i++ )
    nd[i] = INTEGER(r_nd)[i];

  obj = make_mosse_fft(n_fft, nx, dx, nd, flags);

  extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(extPtr, mosse_fft_finalize);

  return extPtr;
}

SEXP r_do_integrate_mosse(SEXP extPtr, SEXP vars, SEXP lambda, SEXP mu, 
		    SEXP drift, SEXP diffusion, SEXP Q, SEXP nt, SEXP dt, 
		    SEXP padding) {
  mosse_fft *obj = (mosse_fft*)R_ExternalPtrAddr(extPtr);
  SEXP ret;
  int nkl = INTEGER(padding)[0], nkr = INTEGER(padding)[1];
  int ndat = LENGTH(lambda);
  double c_dt = REAL(dt)[0];
  int c_nt = INTEGER(nt)[0];
  double *c_lambda=REAL(lambda), *c_mu=REAL(mu);
  double c_drift=REAL(drift)[0], c_diffusion=REAL(diffusion)[0];
  int i, idx, nd;
  if ( obj == NULL )
    error("Corrupt MoSSE integrator: ptr is NULL (are you using multicore?)");
  nd = LENGTH(vars) / obj->nx;
  
  idx = lookup(nd, obj->nd, obj->n_fft);
  if ( idx < 0 )
    error("Failed to find nd = %d\n", nd);

  qf_copy_x_mosse(obj, REAL(vars), nd, 1);

  obj->lambda = REAL(lambda);
  obj->mu = REAL(mu);
  obj->Q = REAL(Q);
    
  for ( i = 0; i < ndat; i++ )
    obj->z[i] = exp(c_dt * (c_lambda[i] - c_mu[i]));
    
  qf_setup_kern_mosse(obj, c_drift, c_diffusion, c_dt, nkl, nkr);

  do_integrate_mosse(obj, c_nt, idx);

  obj->lambda = NULL;
  obj->mu = NULL;

  PROTECT(ret = allocMatrix(REALSXP, obj->nx, nd));
  qf_copy_x_mosse(obj, REAL(ret), nd, 0);
  UNPROTECT(1);
  
  return ret;
}

/* This does the memory allocation and plans the FFT transforms */
mosse_fft* make_mosse_fft(int n_fft, int nx, double dx, int *nd, 
			    int flags) {
  mosse_fft *obj = calloc(1, sizeof(mosse_fft));
  int ny = (((int)floor(nx/2)) + 1);
  int i, max_nd=1;
  for ( i = 0; i < n_fft; i++ )
    if ( nd[i] > max_nd )
      max_nd = nd[i];

  obj->n_fft = n_fft;
  obj->nx = nx;
  obj->ny = ny;
  obj->dx = dx;
  obj->nd = nd;

  obj->x   = fftw_malloc(max_nd *  nx    * sizeof(double));
  obj->y   = fftw_malloc(max_nd * (ny+1) * sizeof(fftw_complex));

  obj->z   = (double*)calloc(nx, sizeof(double));
  obj->wrk = (double*)calloc(nx, sizeof(double));
  obj->wrkd = (double*)calloc(max_nd *  nx, sizeof(double));

  obj->fft = (rfftw_plan_real**)calloc(n_fft, sizeof(rfftw_plan_real*));

  for ( i = 0; i < n_fft; i++ ) {
    obj->fft[i] = make_rfftw_plan_real(nd[i], nx, DIR_COLS,
				       obj->x, obj->y, flags);
  }
  
  /* Brownian kernel */
  obj->kern_x = fftw_malloc(nx     * sizeof(double));
  obj->kern_y = fftw_malloc((ny+1) * sizeof(fftw_complex));
  obj->kernel = make_rfftw_plan_real(1, nx, DIR_COLS, 
				     obj->kern_x, obj->kern_y, flags);
  
  return obj;
}

static void mosse_fft_finalize(SEXP extPtr) {
  mosse_fft *obj = (mosse_fft*)R_ExternalPtrAddr(extPtr);
  int i;
  /* Rprintf("Cleaning up\n"); */

  for ( i = 0; i < obj->n_fft; i++ ) {
    fftw_destroy_plan(obj->fft[i]->plan_f);
    fftw_destroy_plan(obj->fft[i]->plan_b);
  }
  free(obj->fft);
  free(obj->nd);

  fftw_free(obj->x);
  fftw_free(obj->y);

  free(obj->z);
  free(obj->wrk);
  free(obj->wrkd);

  fftw_destroy_plan(obj->kernel->plan_f);
  fftw_destroy_plan(obj->kernel->plan_b);

  fftw_free(obj->kern_x);
  fftw_free(obj->kern_y);

  free(obj);
}


void qf_copy_x_mosse(mosse_fft *obj, double *x, int nd, int copy_in) {
  int i, n = obj->nx * nd;
  double *fft_x = obj->x;
  if ( copy_in )
    for ( i = 0; i < n; i++ )
      fft_x[i] = x[i];
  else {
    for ( i = 0; i < n; i++ ) {
      x[i] = fft_x[i];
    }
  }
}

void qf_setup_kern_mosse(mosse_fft *obj, double drift, double diffusion,
		   double dt, int nkl, int nkr) {
  const int nx = obj->nx;  
  int i;
  double x, *kern_x=obj->kern_x, tot=0, dx=obj->dx;
  double mean, sd;
  
  obj->nkl  = nkl;
  obj->nkr  = nkr;
  obj->npad = nkl + 1 + nkr;
  obj->ndat = nx - obj->npad;
  obj->drift = drift;
  obj->diffusion = diffusion;

  tot   = 0;
  mean  = - dt * drift;
  sd = sqrt(dt * diffusion);

  for ( i = 0, x = 0; i <= nkr; i++, x += dx )
    tot += kern_x[i] = dnorm(x, mean, sd, 0)*dx;
  for ( i = nkr + 1; i < nx - nkl; i++ )
    kern_x[i] = 0;
  for ( i = nx - nkl, x = -nkl*dx; i < nx; i++, x += dx )
    tot += kern_x[i] = dnorm(x, mean, sd, 0)*dx;

  for ( i = 0;        i <= nkr; i++ ) kern_x[i] /= tot;
  for ( i = nx - nkl; i < nx;   i++ ) kern_x[i] /= tot;

  fftw_execute(obj->kernel->plan_f);
}

void do_integrate_mosse(mosse_fft *obj, int nt, int idx) {
  int i, nkl=obj->nkl;
  for ( i = 0; i < nt; i++ ) {
    propagate_t_mosse(obj, idx);
    propagate_x_mosse(obj, idx);
    if ( ISNAN(obj->x[nkl]) )
      error("Integration failure at step %d\n", i);
  }
}

/* Lower level functions */
void propagate_t_mosse(mosse_fft *obj, int idx) {
    int ix, id, ik, nx=obj->nx, ndat=obj->ndat, nd=obj->nd[idx], nk=nd-1, nk2=nk*nk;
  double *vars=obj->x, *d, *dd = obj->wrk;
  double e, tmp1, tmp2, lambda_x, mu_x, z_x, Q_x, d_x;

  for ( ix = 0; ix < ndat; ix++ ) {
    lambda_x = obj->lambda[ix];
    mu_x     = obj->mu[ix];
    z_x      = obj->z[ix];
    e        = vars[ix];
    
    /* Update the E values */
    tmp1 = mu_x - lambda_x * e;
    tmp2 = z_x * (e - 1);
    vars[ix] = (tmp1 + tmp2 * mu_x) / (tmp1 + tmp2 * lambda_x);

    tmp1 = (lambda_x - mu_x) / 
      (z_x*lambda_x - mu_x + (1 - z_x)*lambda_x*e);
    /* Here is the D scaling factor */
    dd[ix] = z_x * tmp1 * tmp1;
  }

  /* Update the D values */
  for ( id = 1; id < nd; id++ ) {
    d = obj->x + nx * id;

    for ( ix = 0; ix < ndat; ix++ ) {
      if ( d[ix] < 0 )
	d[ix] = 0;
      else
	d[ix] *= dd[ix];
	}
  }
  
  for ( id = 1; id < nd; id++ ) {
  	d = obj->wrkd + nx * id;
  	
  	for ( ix = 0; ix < ndat; ix++ ) {
  		d[ix] = 0;
  		
  		for ( ik = 0; ik < nk; ik++) {
  			Q_x = obj->Q[ ix * nk2 + nk * (id-1) + ik ];
  			d_x = obj->x[ nx * (ik+1) + ix ];
  			d[ix] += d_x * Q_x;
  		}
  	}
  }
      
      for ( id = 1; id < nd; id++ ) {
        for ( ix = 0; ix < ndat; ix++ )
             obj->x[nx * id + ix] = obj->wrkd[nx * id + ix];
      }
  }

void propagate_x_mosse(mosse_fft *obj, int idx) {
  double *x, *dd, *wrk = obj->wrk;
  int i, id, nx = obj->nx;
  int nkl = obj->nkl, nkr = obj->nkr, npad = obj->npad;
  int nd=obj->nd[idx];
  
  x = obj->x + 0;
  for ( i = 0; i < nkl; i++ )
    wrk[i] = x[i];
  for ( i = nx-npad-nkr; i < nx - npad; i++ )
    wrk[i] = x[i];
  
  convolve(obj->fft[idx], obj->kern_y);
  
  x = obj->x + 0;
  for ( i = 0; i < nkl; i++ )
    x[i] = wrk[i];
  for ( i = nx-npad-nkr; i < nx - npad; i++ )
    x[i] = wrk[i];
  for ( i = nx - npad; i < nx; i++ ) 
      x[i] = 0; 
    
  for ( id = 1; id < nd; id++ ) {
    x = obj->x + (obj->nx)*id;
    dd = obj->wrkd + (obj->nx)*id;
    for ( i = 0; i < nkl; i++ ) 
      x[i] = dd[i];
    for ( i = nx-npad-nkr; i < nx - npad; i++ )
      x[i] = dd[i];  
    for ( i = nx - npad; i < nx; i++ ) 
      x[i] = 0;  
  }
  
}

#else
SEXP r_do_integrate_mosse(SEXP extPtr, SEXP vars, SEXP star, SEXP b, SEXP mu, SEXP lambda, 
		    SEXP drift, SEXP diffusion, SEXP nt, SEXP dt, 
		    SEXP padding) {
  Rf_error("FFTW support not included");
}
SEXP r_make_mosse_fft(SEXP r_nx, SEXP r_dx, SEXP r_nd, SEXP r_flags) {
  Rf_error("FFTW support not included");
}
#endif
