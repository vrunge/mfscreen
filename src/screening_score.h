#ifndef MFSCREEN_SCREENING_SCORE_H
#define MFSCREEN_SCREENING_SCORE_H

#include <Rcpp.h>

// Squared reciprocal model-free screening statistic for one predictor.
// This is the C++ analogue of screening_scoreR() in the supplied R code.
double screening_score(SEXP x, SEXP y, double tol);

#endif
