use GSL.Filter;
use GSL.RanDist;

const N = 1000; // length of time series
const K:uint(64) = 7; // window size
const f = 5.0; // Frequency of square wave

// Vectors for time series -- these use GSL vectors and so shall we.
var t = new owned GSLVector(N),
  x = new owned GSLVector(N),
  y_median = new owned GSLVector(N),
  y_rmedian = new owned GSLVector(N);

// Generate the input signal
{
  var r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,0);
  for ii in 0.. #N {
    var ti = ii:real/(N-1.0);
    var tmp = sin(2.0*M_PI*f*ti);
    var xi = if (tmp >= 0.0) then 1.0 else -1.0;
    var ei = gsl_ran_gaussian(r, 0.1);

    t[ii] = ti;
    x[ii] = xi + ei;
  }
  gsl_rng_free(r);
}

// Now filter
var median_p = gsl_filter_median_alloc(K);
var rmedian_p = gsl_filter_rmedian_alloc(K);

gsl_filter_median(GSL_FILTER_END_PADVALUE:uint(32),
                  ~x, ~y_median, median_p);
gsl_filter_rmedian(GSL_FILTER_END_PADVALUE:uint(32),
                  ~x, ~y_rmedian, rmedian_p);

for ii in 0.. #N do 
  writef("%dr %dr %dr %dr\n",
         t[ii],
         x[ii],
         y_median[ii],
         y_rmedian[ii]);

gsl_filter_median_free(median_p);
gsl_filter_rmedian_free(rmedian_p);