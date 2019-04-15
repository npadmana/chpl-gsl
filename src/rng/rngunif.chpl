use GSL.RNG;

config const RNGType="mt19937";
config const iseed : uint(64) = 0;

var r : c_ptr(gsl_rng);
if RNGType=="mrg" {
  r = gsl_rng_alloc(gsl_rng_mrg);
} else {
  r = gsl_rng_alloc(gsl_rng_mt19937);
}
gsl_rng_set(r, iseed);

for i in 1..10 do writef("%.5dr\n", gsl_rng_uniform(r));

gsl_rng_free(r);



