use GSL.Sorting;
use GSL.RNG;

const N = 100000;
const K = 5;

var x : [0.. #N] real;
var small : [0.. #K] real;

var r = gsl_rng_alloc(gsl_rng_mt19937);
gsl_rng_set(r,0);
// Do this serially to match output
for xi in x do xi = gsl_rng_uniform(r);
gsl_rng_free(r);

gsl_sort_smallest(~small, K:uint(64), ~x, 1, N:uint(64));

writef("%u smallest values from %u\n",K,N);
for ii in small.domain {
  writef("%u: %.18dr\n",ii, small[ii]);
}


