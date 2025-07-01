
#ifndef INFERNO_MCMC_H
#define INFERNO_MCMC_H

#include <R.h>
#include <Rinternals.h>

SEXP metrop(SEXP func1, SEXP initial, SEXP nbatch, SEXP blen, SEXP nspac,
    SEXP scale, SEXP func2, SEXP debug, SEXP rho1, SEXP rho2);

SEXP temper(SEXP func1, SEXP initial, SEXP neighbors, SEXP nbatch,
    SEXP blen, SEXP nspac, SEXP scale, SEXP func2, SEXP debug,
    SEXP parallel, SEXP rho1, SEXP rho2);

SEXP initseq(SEXP x);

void olbm(double *x, int *nin, int *pin, int *lin, double *mean,
    double *var, int *nocalcin);

#endif /* INFERNO_MCMC_H */

