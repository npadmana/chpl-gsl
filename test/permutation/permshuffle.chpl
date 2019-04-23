use GSL.Permutations;
use GSL.RanDist;

const N : size_t = 10;

var r = gsl_rng_alloc(gsl_rng_mt19937);
gsl_rng_set(r, 0);

var p = gsl_permutation_alloc(N);
var q = gsl_permutation_alloc(N);

writef("initial permutation:");
gsl_permutation_init(p);
print_perm(p);

// Hardcode size_t
gsl_ran_shuffle(r, gsl_permutation_data(p), N, 8);
writef(" random permutation:");
print_perm(p);

gsl_permutation_inverse(q,p);
writef("inverse permutation:");
print_perm(q);

gsl_permutation_free(p);
gsl_permutation_free(q);
gsl_rng_free(r);


proc print_perm(p) {
  var data = gsl_permutation_data(p);
  var n = gsl_permutation_size(p);
  for i in 0.. #n do writef(" %i",data[i]);
  writef("\n");
}