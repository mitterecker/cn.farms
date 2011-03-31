#include <Rdefines.h>

/* laplace.cpp */
SEXP momentsGauss(SEXP it, SEXP eps1S, SEXP eps2S, SEXP aS, SEXP bS, SEXP cS, SEXP sigmaZS, SEXP normfactS, SEXP evalProb, SEXP methodS);

/* sparse_farms.c */
SEXP sparseFarmsC(SEXP xS, SEXP cycS, SEXP XXS, SEXP nnS);
SEXP getL(SEXP pointer);
SEXP getEss(SEXP pointer);
SEXP getLap(SEXP pointer);
SEXP deinit(SEXP pointer);

