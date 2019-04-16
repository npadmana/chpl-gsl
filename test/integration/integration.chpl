use GSL.Integration;

proc main() {
  var w = gsl_integration_workspace_alloc(1000);
  defer gsl_integration_workspace_free(w);

  const expected = -4.0;
  var alpha : real = 1.0;

  var F = mkGSLFunc(c_ptrTo(func1), alpha);

  var result, error : c_double;
  gsl_integration_qags(~F, 0.0, 1.0, 0.0, 1.0e-7,
                       1000, w, ~result, ~error); 

  writef("result          = % .18dr\n", result);
  writef("exact result    = % .18dr\n", expected);
  writef("estimated error = % .18dr\n", error);
  writef("actual error    = % .18dr\n", result - expected);
  writef("intervals       = %u\n", (w.deref()).size);
}



proc func1(x : c_double, params : c_void_ptr) : c_double {
  const alpha = params.deref(c_double);
  const ret = log(alpha * x)/sqrt(x);
  return ret;
}
