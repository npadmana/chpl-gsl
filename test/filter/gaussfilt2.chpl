use GSL.Filter;
use GSL.RanDist;

const N = 1000; // length of time series
const K:uint(64) = 61; // window size
const alpha = 3.0;

// Vectors for time series -- these use GSL vectors and so shall we.
var x = new owned GSLVector(N),
  y = new owned GSLVector(N),
  dy = new owned GSLVector(N),
  d2y = new owned GSLVector(N);

// Generate the input signal
{
  var r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,0);
  var sum = 0.0;
  for ii in 0.. #N {
    var xi = if ii > N/2 then 0.5 else 0.0;
    var ei = gsl_ran_gaussian(r, 0.1);
    x[ii] = xi+ei;
  }
  gsl_rng_free(r);
}

// Now filter
var gauss_p = gsl_filter_gaussian_alloc(K);
gsl_filter_gaussian(GSL_FILTER_END_PADVALUE:uint(32),
                    alpha, 0, ~x, ~y, gauss_p);
gsl_filter_gaussian(GSL_FILTER_END_PADVALUE:uint(32),
                    alpha, 1, ~x, ~dy, gauss_p);
gsl_filter_gaussian(GSL_FILTER_END_PADVALUE:uint(32),
                    alpha, 2, ~x, ~d2y, gauss_p);

for ii in 0.. #N {
  var dxi : real;
  var xi = x[ii];
  if ii == 0 {
    dxi = x[1]-xi;
  } else if ii == N-1 {
    dxi = x[ii]-x[ii-1];
  } else {
    dxi = (x[ii+1] - x[ii-1])*0.5;
  }

  writef("%.12er %.12er %.12er %.12er %.12er\n",
         xi,
         y[ii],
         dxi,
         dy[ii],
         d2y[ii]);
}
gsl_filter_gaussian_free(gauss_p);