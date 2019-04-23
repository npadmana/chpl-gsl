use GSL.Combinations;

const n : size_t = 4;

writef("All subsets of {0,1,2,3} by size:\n") ;
for i in 0..4 {
  var c = gsl_combination_calloc(4, i : uint(64));
  do {
    writef("{");
    var data = gsl_combination_data(c);
    for j in 0.. #i do writef(" %i",data[j]);
    writef(" }\n");
  } while (gsl_combination_next(c) == GSL_SUCCESS);
  gsl_combination_free(c);
}
