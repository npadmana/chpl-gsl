module demo {
  use GSL;

  record quadratic_params {
    var a, b, c: real;
  }

  proc quadratic(x : real, params : c_void_ptr) : real {
    ref p = params.deref(quadratic_params);

    return (p.a * x + p.b)*x + p.c;
  }

  proc quadratic_deriv(x : real, params : c_void_ptr) : real {
    ref p = params.deref(quadratic_params);

    return 2.0*p.a*x + p.b;
  }

  proc quadratic_fdf(x : real, params : c_void_ptr,
                     y : c_ptr(c_double),
                     dy : c_ptr(c_double)) {
    ref p = params.deref(quadratic_params);

    y.deref()=  (p.a * x + p.b)*x + p.c;
    dy.deref()= 2.0*p.a*x + p.b;
  }

}