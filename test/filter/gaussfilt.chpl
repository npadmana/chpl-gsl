use GSL.Filter;
use GSL.RanDist;

const N = 500; // length of time series
const K:uint(64) = 51; // window size
const alpha = [0.5, 3.0, 10.0];

// Vectors for time series -- these use GSL vectors and so shall we.
var x = new owned GSLVector(N),
  y1 = new owned GSLVector(N),
  y2 = new owned GSLVector(N),
  y3 = new owned GSLVector(N);

// Vectors for kernels
var k1 = new owned GSLVector(K),
  k2 = new owned GSLVector(K),
  k3 = new owned GSLVector(K);

// Generate the input signal
{
  var r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,0);
  var sum = 0.0;
  for ii in 0.. #N {
    var ui = gsl_ran_gaussian(r,1.0);
    sum += ui;
    x[ii] = sum;
  }
  gsl_rng_free(r);
}

// Compute the kernels without normalization
const order : size_t = 0;
const normalize : c_int = 0;
gsl_filter_gaussian_kernel(alpha[1], order, normalize, ~k1);
gsl_filter_gaussian_kernel(alpha[2], order, normalize, ~k2);
gsl_filter_gaussian_kernel(alpha[3], order, normalize, ~k3);

// Print kernels
for ii in 0.. #K do writef("%er %er %er\n",k1[ii],k2[ii],k3[ii]);
writef("\n\n");

// Now filter
var gauss_p = gsl_filter_gaussian_alloc(K);
gsl_filter_gaussian(GSL_FILTER_END_PADVALUE:uint(32),
                    alpha[1], 0, ~x, ~y1, gauss_p);
gsl_filter_gaussian(GSL_FILTER_END_PADVALUE:uint(32),
                    alpha[2], 0, ~x, ~y2, gauss_p);
gsl_filter_gaussian(GSL_FILTER_END_PADVALUE:uint(32),
                    alpha[3], 0, ~x, ~y3, gauss_p);

for ii in 0.. #N do writef("%.12er %.12er %.12er %.12er\n",x[ii],y1[ii],y2[ii],y3[ii]);
gsl_filter_gaussian_free(gauss_p);