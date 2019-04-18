// BSpline example program.
// In the original output, the chi2 values got written to stderr and
// were not in the file. I've added these to the file; the original
// data are still here as well.

// Unlike the original version, this uses GSL matrices and vectors instead
// of Chapel matrices and vectors. This might significantly simplify some
// code, which no longer has to dance around views.

use GSL.BSpline;
use GSL.RanDist;
use GSL.LinearFit;
use GSL.Statistics;

config const N:size_t=200; // Number of data points
config const NCOEFFS:size_t=12; // Number of fit coefficients
config const ORDER:size_t=4;
const NBREAK:size_t=NCOEFFS + 2 - ORDER;

// Fill the data vectors
var x = new owned GSLVector(N),
  y = new owned GSLVector(N),
  w = new owned GSLVector(N);
var r = gsl_rng_alloc(gsl_rng_mt19937);
gsl_rng_set(r,0);
for ii in 0.. #N {
  x[ii] = (15.0/(N-1))*ii;
  y[ii] = cos(x[ii])*exp(-0.1*x[ii]);
  const sigma = 0.1*y[ii];
  y[ii] += gsl_ran_gaussian(r, sigma);
  w[ii] = 1.0/(sigma**2);
  writef("%.6dr %.6dr\n",x[ii],y[ii]);
}

var X = new owned GSLMatrix(N, NCOEFFS);
var B = new owned GSLVector(NCOEFFS);
var bw = gsl_bspline_alloc(4, NBREAK);

/* use uniform breakpoints of [0, 15] */
gsl_bspline_knots_uniform(0.0, 15.0, bw);

// Construct the fit matrix
for ii in 0..#N {
  gsl_bspline_eval(x[ii], B.p, bw);
  forall jj in 0.. #NCOEFFS do X[ii,jj] = B[jj];
}
  
// Now do the fit
var c = new owned GSLVector(NCOEFFS);
var cov = new owned GSLMatrix(NCOEFFS, NCOEFFS);
var mw = gsl_multifit_linear_alloc(N, NCOEFFS);

// Get GSL views of matrices
var chisq:real;
gsl_multifit_wlinear(~X, ~w, ~y, ~c, ~cov,
                     ~chisq,
                     mw);

const dof = N - NCOEFFS;
const tss = gsl_stats_wtss(w.pdata, 1, y.pdata, 1, N);
const Rsq = 1.0 - chisq/tss;
writef("chisq/dof = %er, Rsq = %dr\n\n\n",chisq/dof, Rsq);

var xi, yi, yerr:real;
xi = 0.0;
while (xi < 15.0) {
  gsl_bspline_eval(xi, B.p, bw);
  gsl_multifit_linear_est(~B, ~c, ~cov, ~yi, ~yerr);
  writef("%.6dr %.6dr\n",xi,yi);
  xi += 0.1;
}

gsl_rng_free(r);
gsl_bspline_free(bw);
gsl_multifit_linear_free(mw);

