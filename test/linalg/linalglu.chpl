use GSL.LinearAlgebra;

var a_data = reshape([ 0.18, 0.60, 0.57, 0.96,
                       0.41, 0.24, 0.99, 0.58,
                       0.14, 0.30, 0.97, 0.66,
                       0.51, 0.13, 0.19, 0.85 ],
                     {0.. #4, 0.. #4});
var b_data = [ 1.0, 2.0, 3.0, 4.0 ];

var m = GSLView(a_data);
var b = GSLView(b_data);

var x = new owned GSLVector(4);

var s : c_int;

var p = gsl_permutation_alloc(4);

// Integers need an explicit cast
gsl_linalg_LU_decomp(~m, p, c_ptrTo(s));
gsl_linalg_LU_solve(~m, p, ~b, ~x);

writef("x = \n");
for ii in x.D do writef("%r\n",x[ii]);

gsl_permutation_free(p);