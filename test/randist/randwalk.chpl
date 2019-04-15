use GSL.RanDist;

config const n=10;
config const seed=123;

var r = gsl_rng_alloc(gsl_rng_mt19937);
gsl_rng_set(r, 123);

var x = 0.0, y=0.0;
writef("%8.5dr %8.5dr\n",x,y);
var dx, dy : c_double;

for i in 0.. #n {
  gsl_ran_dir_2d(r, ~dx, ~dy);
  x += dx;
  y += dy;
  writef("%8.5dr %8.5dr\n",x,y);
}
gsl_rng_free(r);