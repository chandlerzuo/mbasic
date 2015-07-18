#ifndef _csnet_RCPP_HELLO_WORLD_H
#define _csnet_RCPP_HELLO_WORLD_H

#include <Rcpp.h>
//#include <stdlib.h>
#include <time.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
using namespace Rcpp;
RcppExport SEXP mcmc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) ;

RcppExport SEXP madbayes(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP madbayes_theta(SEXP, SEXP,SEXP,SEXP,SEXP, SEXP, SEXP);

double ComputeLoss(IntegerVector, IntegerMatrix, NumericMatrix, NumericMatrix, NumericMatrix, NumericMatrix, IntegerVector, IntegerVector, double, double);
double ComputeLoss_theta(IntegerVector, IntegerMatrix, NumericMatrix, NumericMatrix, NumericMatrix);

RcppExport SEXP RSolveW(SEXP, SEXP, SEXP, SEXP);
NumericVector SolveW(NumericVector, double, double, double);

RcppExport SEXP madbayes_init(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP madbayes_init_kmeanspp(SEXP, SEXP, SEXP);
int SampleCentroid(NumericVector);

#endif
