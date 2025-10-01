/*
 * These are the protracted_representative equations, implemented in c
 */
#include <R.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_protractedreptip(double *pars, const double *y, double *ydot) {
  double Ei = y[0], Er = y[1], Dr = y[2], Di = y[3], Dir = y[4];
  double b = pars[0], mu = pars[1], lambda = pars[2];

  ydot[0] = -(b+mu) * Ei + b * Ei * Ei + mu;
  ydot[1] = -(b+mu+lambda) * Er + b * Er * Er + mu + lambda * Ei;
  ydot[2] = -(b+mu) * Dr + 2 * b * Dr * Ei + lambda * (Er - Ei);
  ydot[3] = -(b+mu+lambda) * Di + 2 * b * Di * Er + lambda * (Er - Ei + Dr);
  ydot[4] = -(b+mu+lambda) * Dir + 2 * b * Dir * Ei + lambda * (Er - Ei + Dr);
}

/* 
   Wrap these up for GslOde, which allows time and neq to be here too.
   
   TODO: Eventually move the contents of do_derivs_* into these
   functions, as no other backend possible now.
 */
void derivs_protractedreptip_gslode(int neqs, double t, double *pars, 
			 const double *y, double *dydt) {
  do_derivs_protractedreptip(pars, y, dydt);
}
