use GSL.RanDist;

config const n=10;
config const mu=3.0;

var r = gsl_rng_alloc(gsl_rng_mt19937);
gsl_rng_set(r, 123);

for i in 0.. #n {
  const k = gsl_ran_poisson(r, mu);
  writef(" %i", k);
}
writef("\n");
gsl_rng_free(r);