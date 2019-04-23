use GSL.Permutations;

var p = gsl_permutation_alloc(3);

gsl_permutation_init(p);

do {
  var data = gsl_permutation_data(p);
  for i in 0.. #3 do writef(" %i",data[i]);
  writef("\n");
 } while (gsl_permutation_next(p) == GSL_SUCCESS);

gsl_permutation_free(p);