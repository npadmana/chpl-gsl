use GSL.QRNG;

var q = gsl_qrng_alloc(gsl_qrng_sobol, 2);
var v : [0..1]real;

for i in 0.. #1024  {
  gsl_qrng_get(q, ~v);
  writef("%.5dr %.5dr\n", v[0], v[1]);
}

gsl_qrng_free(q);