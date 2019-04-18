use GSL.Chebyshev;

const n=10000;

var cs = gsl_cheb_alloc(40);
var F = mkGSLFunc(c_ptrTo(step));

gsl_cheb_init(cs, ~F, 0.0, 1.0);
for ii in 0.. #n {
  const x= ii/n:real;
  const r10 = gsl_cheb_eval_n(cs,10,x);
  const r40 = gsl_cheb_eval(cs,x);
  writef("%.4dr %.2dr %.6dr %.6dr\n",
         x, step(x, nil : c_void_ptr), r10, r40);
}

gsl_cheb_free(cs);



proc step(x : real, p : c_void_ptr) : real {
  if (x < 0.5) {
    return 0.25;
  } else {
    return 0.75;
  }
}
