use GSL.Polynomials;

var coeff : [0.. #6]real = [-1,0,0,0,0,1];
var z : [0.. #5]complex;

var w = gsl_poly_complex_workspace_alloc(6);
gsl_poly_complex_solve(~coeff, 6, w, (~z):gsl_complex_packed_ptr);
gsl_poly_complex_workspace_free(w);

for ii in z.domain do writef("z%i = %+.18dr %+.18dr\n",ii,z[ii].re,z[ii].im);