use GSL.Deriv;

var F = mkGSLFunc(c_ptrTo(simple));
var result, abserr : real;


writef("f(x) = x^(3/2)\n");

gsl_deriv_central(~F, 2.0, 1e-8, ~result, ~abserr);
writef("x = 2.0\n");
writef("f'(x) = %.10dr +/- %.10dr\n", result, abserr);
writef("exact = %.10dr\n\n", 1.5 * sqrt(2.0));

gsl_deriv_forward(~F, 0.0, 1e-8, ~result, ~abserr);
writef("x = 0.0\n");
writef("f'(x) = %.10dr +/- %.10dr\n", result, abserr);
writef("exact = %.10dr\n", 0.0);


proc simple(x : real, params : c_void_ptr) : real {
  return x**1.5;
}