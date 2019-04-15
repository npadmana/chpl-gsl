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

writef("generator type: %s\n", gsl_rng_name(r):string);
writef("seed = %i\n", iseed);
writef("first value = %i\n", gsl_rng_get(r));

gsl_rng_free(r);



