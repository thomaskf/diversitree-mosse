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
#include "protracted-mutation-eqs-fftC.h"

static void protracted_mutation_fft_finalize(SEXP extPtr);

SEXP r_make_protracted_mutation_fft(SEXP r_nx, SEXP r_dx, SEXP r_nd, SEXP r_flags) {
  protracted_mutation_fft *obj;
  SEXP extPtr;
  int nx = INTEGER(r_nx)[0];
  double dx = REAL(r_dx)[0];
  int n_fft = LENGTH(r_nd);
  int i;
  int flags;
  int *nd = (int*)calloc(n_fft, sizeof(int));
  for ( i = 0; i < n_fft; i++ )
    nd[i] = INTEGER(r_nd)[i];
  
  if ( INTEGER(r_flags)[0] == -1 )
    flags = FFTW_ESTIMATE;
  else if ( INTEGER(r_flags)[0] == 1 )
    flags = FFTW_PATIENT;
  else if ( INTEGER(r_flags)[0] == 2 )
    flags = FFTW_EXHAUSTIVE;
  else
    flags = FFTW_MEASURE;

  obj = make_protracted_mutation_fft(n_fft, nx, dx, nd, flags);

  extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(extPtr, protracted_mutation_fft_finalize);

  return extPtr;
}

SEXP r_do_integrate_protracted_mutation(SEXP extPtr, SEXP vars, SEXP star, SEXP b, SEXP mu, SEXP lambda, 
		    SEXP drift, SEXP diffusion, SEXP nt, SEXP dt, 
		    SEXP padding) {
  protracted_mutation_fft *obj = (protracted_mutation_fft*)R_ExternalPtrAddr(extPtr);
  SEXP ret;
  int nkl = INTEGER(padding)[0], nkr = INTEGER(padding)[1];
  int ndat = LENGTH(lambda);
  int c_star = INTEGER(star)[0];
  double c_dt = REAL(dt)[0], c_nt = INTEGER(nt)[0];
  double *c_b=REAL(b), *c_mu=REAL(mu), *c_lambda=REAL(lambda);
  double c_drift=REAL(drift)[0], c_diffusion=REAL(diffusion)[0];
  int i, idx, nd;
  if ( obj == NULL )
    error("Corrupt ProSSE integrator: ptr is NULL (are you using multicore?)");
  nd = LENGTH(vars) / obj->nx;
  
  idx = lookup(nd, obj->nd, obj->n_fft);
  if ( idx < 0 )
    error("Failed to find nd = %d\n", nd);

  qf_copy_x_protracted_mutation(obj, REAL(vars), nd, 1);

  obj->b = REAL(b);
  obj->mu = REAL(mu);
  obj->lambda = REAL(lambda);
  for ( i = 0; i < ndat; i++ )
    obj->z1[i] = exp(c_dt * (c_b[i] - c_mu[i]));
    obj->z2[i] = exp(c_dt * (- c_lambda[i]));
    
  qf_setup_kern_protracted_mutation(obj, c_drift, c_diffusion, c_dt, nkl, nkr);

  do_integrate_protracted_mutation(obj, c_star, c_nt, idx);

  obj->b = NULL;
  obj->mu = NULL;
  obj->lambda = NULL;

  PROTECT(ret = allocMatrix(REALSXP, obj->nx, nd));
  qf_copy_x_protracted_mutation(obj, REAL(ret), nd, 0);
  UNPROTECT(1);
  
  return ret;
}

/* This does the memory allocation and plans the FFT transforms */
protracted_mutation_fft* make_protracted_mutation_fft(int n_fft, int nx, double dx, int *nd, 
			    int flags) {
  protracted_mutation_fft *obj = calloc(1, sizeof(protracted_mutation_fft));
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

  obj->z1   = (double*)calloc(nx, sizeof(double));
  obj->z2   = (double*)calloc(nx, sizeof(double));
  obj->wrk = (double*)calloc(nx, sizeof(double));

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

static void protracted_mutation_fft_finalize(SEXP extPtr) {
  protracted_mutation_fft *obj = (protracted_mutation_fft*)R_ExternalPtrAddr(extPtr);
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

  free(obj->z1);
  free(obj->z2);
  free(obj->wrk);

  fftw_destroy_plan(obj->kernel->plan_f);
  fftw_destroy_plan(obj->kernel->plan_b);

  fftw_free(obj->kern_x);
  fftw_free(obj->kern_y);

  free(obj);
}


void qf_copy_x_protracted_mutation(protracted_mutation_fft *obj, double *x, int nd, int copy_in) {
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

void qf_copy_ED_protracted_mutation(protracted_mutation_fft *obj, double *x, int idx) {
  int i, nx = obj->nx;
  double *fft_x = obj->x;
  for ( i = 0; i < nx; i++ )
    x[i] = fft_x[i];

  x += nx;
  fft_x = obj->x + idx*nx;

  for ( i = 0; i < nx; i++ )
    x[i] = fft_x[i];
}

void qf_setup_kern_protracted_mutation(protracted_mutation_fft *obj, double drift, double diffusion,
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

void do_integrate_protracted_mutation(protracted_mutation_fft *obj, int star, int nt, int idx) {
  int i, nkl=obj->nkl;
  for ( i = 0; i < nt; i++ ) {
    propagate_t_protracted_mutation(obj, star, idx);
    propagate_x_protracted_mutation(obj, idx);
    if ( ISNAN(obj->x[nkl]) )
      error("Integration failure at step %d\n", i);
  }
}

/* Lower level functions */
void propagate_t_protracted_mutation(protracted_mutation_fft *obj, int star, int idx) {
  int ix, nx=obj->nx, ndat=obj->ndat;
  double *vars=obj->x, *d, *dr, *dd = obj->wrk;
  double e, tmp1, tmp2, b_x, mu_x, lambda_x, z1_x, z2_x;

  for ( ix = 0; ix < ndat; ix++ ) {
    b_x      = obj->b[ix];
    mu_x     = obj->mu[ix];
    lambda_x = obj->lambda[ix];
    z1_x      = obj->z1[ix];
    e        = vars[ix];
    
    /* Update the E values */
    tmp1 = mu_x - b_x * e;
    tmp2 = z1_x * (e - 1);
    vars[ix] = (tmp1 + tmp2 * mu_x) / (tmp1 + tmp2 * b_x);

    tmp1 = (b_x - mu_x) / 
      (z1_x*b_x - mu_x + (1 - z1_x)*b_x*e);
    /* Here is the D scaling factor */
    dd[ix] = z1_x * tmp1 * tmp1;
  }

  /* Update the D values */
  d = obj->x + nx;
  dr = obj->x + 2 * nx;
  if (star == 1) {
    for ( ix = 0; ix < ndat; ix++ ) {
    	z2_x = obj->z2[ix];
    	tmp1 = dd[ix] * z2_x;
      	if ( dr[ix] < 0 ) {
			dr[ix] = 0;
      	} else {
			dr[ix] *= tmp1;
	  	}
	  	if ( d[ix] < 0 ) {
			d[ix] = 0;
      	} else {
			d[ix] *= tmp1;
	  	}
  	}
  } else if (star == 2){
    for ( ix = 0; ix < ndat; ix++ ) {
    	z2_x = obj->z2[ix];
    	tmp1 = dr[ix];
    	if ( dr[ix] < 0 ) {
			dr[ix] = 0;
      	} else {
			dr[ix] *= dd[ix];
	  	}
	  	if ( d[ix] < 0 ) {
			d[ix] = 0;
      	} else {
			d[ix] = dr[ix] + (d[ix] - tmp1) * dd[ix] * z2_x;
	  	}
	  	
  	}
  }
}

void propagate_x_protracted_mutation(protracted_mutation_fft *obj, int idx) {
  double *x = obj->x, *z1 = obj->z1, *z2 = obj->z2, *wrk = obj->wrk;
  int i, id, nx = obj->nx;
  int nkl = obj->nkl, nkr = obj->nkr, npad = obj->npad;
  int nd=obj->nd[idx];

  /* TODO: I am not sure if these are flipped nkl/nkr */
  for ( i = 0; i < nkl; i++ )
    z1[i] = x[i];
  for ( i = nx-npad-nkr; i < nx - npad; i++ )
    z1[i] = x[i];
  for ( i = nx; i < nkl + nx; i++ )
    z2[i] = x[i];
  for ( i = 2 * nx-npad-nkr; i < 2 * nx - npad; i++ )
    z2[i] = x[i];
  for ( i = 2* nx; i < nkl + 2 * nx; i++ )
    wrk[i] = x[i];
  for ( i = 3 * nx-npad-nkr; i < 3 * nx - npad; i++ )
    wrk[i] = x[i];

  convolve(obj->fft[idx], obj->kern_y);

  for ( i = 0; i < nkl; i++ )
    x[i] = z1[i];
  for ( i = nx-npad-nkr; i < nx - npad; i++ )
    x[i] = z1[i];
  for ( i = nx; i < nkl + nx; i++ )
    x[i] = z2[i];
  for ( i = 2 * nx-npad-nkr; i < 2 * nx - npad; i++ )
    x[i] = z2[i];
  for ( i = 2* nx; i < nkl + 2 * nx; i++ )
    x[i] = wrk[i];
  for ( i = 3 * nx-npad-nkr; i < 3 * nx - npad; i++ )
    x[i] = wrk[i];

  /* Zeroing takes a little more work, now.  We might be able to get
     away with just zeroing the E though */
  for ( id = 0; id < nd; id++ ) {
    x = obj->x + (obj->nx)*(id+1) - npad;
    for ( i = 0; i < npad; i++ ) 
      x[i] = 0;
  }
}

#else
SEXP r_do_integrate_protracted_mutation(SEXP extPtr, SEXP vars, SEXP star, SEXP b, SEXP mu, SEXP lambda, 
		    SEXP drift, SEXP diffusion, SEXP nt, SEXP dt, 
		    SEXP padding) {
  Rf_error("FFTW support not included");
}
SEXP r_make_protracted_mutation_fft(SEXP r_nx, SEXP r_dx, SEXP r_nd, SEXP r_flags) {
  Rf_error("FFTW support not included");
}
#endif
