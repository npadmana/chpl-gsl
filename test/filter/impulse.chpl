use GSL.Filter;
use GSL.RanDist;

const N = 1000; // length of time series
const K:uint(64) = 25; // window size
const t = 4.0; // # of scale factors for outlier detection

// Vectors for time series -- these use GSL vectors and so shall we.
var x = new owned GSLVector(N),
  y = new owned GSLVector(N),
  xmedian = new owned GSLVector(N),
  xsigma = new owned GSLVector(N),
  ioutlier = new owned GSLVector(c_int, N);

// Generate the input signal
{
  var r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,0);
  for ii in 0.. #N {
    var ti = ii:real/N;
    var xi = 10.0*sin(2.0*M_PI*ti);
    var ei = gsl_ran_gaussian(r, 2.0);
    var u = gsl_rng_uniform(r);
    var outlier : real = if (u < 0.01) then 15.0*sgn(ei) else 0.0;
    x[ii] = xi+ei+outlier;
  }
  gsl_rng_free(r);
}

// Now impulse detect
var w = gsl_filter_impulse_alloc(K);
var noutlier : uint(64);
gsl_filter_impulse(GSL_FILTER_END_TRUNCATE:uint(32),
                   GSL_FILTER_SCALE_QN:uint(32),
                   t, ~x, ~y, ~xmedian, ~xsigma,
                   c_ptrTo(noutlier), ~ioutlier, w);


for ii in 0.. #N do 
      writef("%u %dr %dr %dr %dr %i\n",
             ii,
             x[ii],
             y[ii],
             xmedian[ii]+t*xsigma[ii],
             xmedian[ii]-t*xsigma[ii],
             ioutlier[ii]);

gsl_filter_impulse_free(w);
