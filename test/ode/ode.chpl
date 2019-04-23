use GSL.ODE;

config const choice="init";

var mu = 10.0;
select choice {
  when "init" do initval();
  when "low" do lowlevel();
  when "fixed" do fixed();
}

proc fixed() {
  // The initval version of the code
  var sys : gsl_odeiv2_system;
  sys.function = c_ptrTo(odefunc);
  sys.jacobian = c_ptrTo(odejac);
  sys.dimension = 2;
  sys.params = ~mu;

  var d = gsl_odeiv2_driver_alloc_y_new(~sys,
                                    gsl_odeiv2_step_rk4,
                                    1.0e-3,
                                    1.0e-8,
                                    1.0e-8);
  defer gsl_odeiv2_driver_free(d);

  var y = [1.0, 0.0];
  var t : real;
  for i in 1..100 {
    var status = gsl_odeiv2_driver_apply_fixed_step(d,
                                                    ~t,
                                                    1.0e-3,
                                                    1000,
                                                    ~y);
    if status != GSL_SUCCESS {
      writef("error, return value=%i\n", status);
      break;
    }
    writef("%.5er %.5er %.5er\n", t, y[1], y[2]);
  }

}

proc lowlevel() {
  var sys : gsl_odeiv2_system;
  sys.function = c_ptrTo(odefunc);
  sys.jacobian = c_ptrTo(odejac);
  sys.dimension = 2;
  sys.params = ~mu;

  var s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd, 2);
  defer gsl_odeiv2_step_free(s);
  var c = gsl_odeiv2_control_y_new(1.0e-6, 0.0);
  defer gsl_odeiv2_control_free(c);
  var e = gsl_odeiv2_evolve_alloc(2);
  defer gsl_odeiv2_evolve_free(e);

  var y = [1.0, 0.0];
  var t : real = 0.0;
  const t1 = 100.0;
  var h : real = 1.0e-6;

  while (t < t1) {
    var status = gsl_odeiv2_evolve_apply(e, c, s,
                                         ~sys,
                                         ~t,
                                         t1,
                                         ~h,
                                         ~y);
    if status != GSL_SUCCESS then break;

    writef("%.5er %.5er %.5er\n", t, y[1], y[2]);
  }
}


proc initval() {
  // The initval version of the code
  var sys : gsl_odeiv2_system;
  sys.function = c_ptrTo(odefunc);
  sys.jacobian = c_ptrTo(odejac);
  sys.dimension = 2;
  sys.params = ~mu;

  var d = gsl_odeiv2_driver_alloc_y_new(~sys,
                                    gsl_odeiv2_step_rk8pd,
                                    1.0e-6,
                                    1.0e-6,
                                    0.0);
  defer gsl_odeiv2_driver_free(d);

  var y = [1.0, 0.0];
  var t : real;
  const t1 = 100.0;
  for i in 1..100 {
    const ti = i*t1/100.0;
    var status = gsl_odeiv2_driver_apply(d, ~t, ti, ~y);
    if status != GSL_SUCCESS {
      writef("error, return value=%i\n", status);
      break;
    }
    writef("%.5er %.5er %.5er\n", t, y[1], y[2]);
  }

}



/* Define functions for function and jacobian */
proc odefunc(t : real, in y : c_ptr(real), f : c_ptr(real),
             params : c_void_ptr) : c_int {
  ref mu = params.deref(real);
  f[0] = y[1];
  f[1] = -y[0] - mu*y[1]*(y[0]**2 - 1);
  return GSL_SUCCESS;
}

proc odejac(t : real, in y : c_ptr(real),
            dfdy : c_ptr(real),
            dfdt : c_ptr(real), params : c_void_ptr) : c_int {
  ref mu = params.deref(real);
  // 2*2 array -- do this by hand.
  // A more general case would do this
  // by defining a matrix
  dfdy[0*2+0] = 0.0;
  dfdy[0*2+1] = 1.0;
  dfdy[1*2+0] = -2.0*mu*y[0]*y[1] - 1.0;
  dfdy[1*2+1] = -mu*(y[0]**2  - 1.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

