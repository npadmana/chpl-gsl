use GSL.SpecFun;
use GSL.Integration;

config var m : c_int = 10;
config const n : size_t = 6;

proc main() {
  var expected : real;
  if (m % 2 == 0) {
    expected = M_SQRTPI + gsl_sf_gamma(0.5*(1+m));
  } else {
    expected = M_SQRTPI;
  }

  var F = mkGSLFunc(c_ptrTo(func1), m);
  //var F : gsl_function;
  //F.function = c_ptrTo(func1);
  //F.params = c_ptrTo(m):c_void_ptr;

  var w = gsl_integration_fixed_alloc(
                                      gsl_integration_fixed_hermite,
                                      n, 0.0, 1.0, 0.0, 0.0);
  defer gsl_integration_fixed_free(w);

  var result : c_double;
  gsl_integration_fixed(~F, ~result, w);

  writef("m               = %u\n", m);
  writef("intervals       = %u\n", gsl_integration_fixed_n(w));
  writef("result          = % .18dr\n", result);
  writef("exact result    = % .18dr\n", expected);
  writef("actual error    = % .18dr\n", result - expected);
}



proc func1(x : c_double, params : c_void_ptr) : c_double {
  const m1 = params.deref(c_int);
  return x**m1 + 1.0;
}
