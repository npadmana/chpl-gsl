/* Example using derivatives for minimization.

   This version uses conjugate gradients. */

use GSL.Minimize;

// Starting point
var x : [0.. #2] real = [5.0, 7.0];

// Parameters of parabaloid
var pp = new ParamStruct(1.0, 2.0, 10.0, 20.0, 30.0);

// Set up user function and minimizer
var minex_func : gsl_multimin_function_fdf;
minex_func.n = 2;
minex_func.f = c_ptrTo(my_f);
minex_func.df = c_ptrTo(my_df);
minex_func.fdf = c_ptrTo(my_fdf);
minex_func.params = c_ptrTo(pp);

var s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, 2);
// Get GSL vector views on arrays
var xview = GSLView(x);
gsl_multimin_fdfminimizer_set(s, ~minex_func,
                              ~xview, 0.01, 1.0e-4);

// Iterations
var status : c_int;
var itnum = 0;

do {
  itnum += 1;
  status = gsl_multimin_fdfminimizer_iterate(s);

  if status then break;

  var grad = gsl_multimin_fdfminimizer_gradient(s);
  status = gsl_multimin_test_gradient(grad, 1.0e-3);

  if status==GSL_SUCCESS then writef("Minimum found at:\n");
  var xminp = gsl_multimin_fdfminimizer_x(s);
  var fmin = gsl_multimin_fdfminimizer_minimum(s);
  writef("%5i %.5dr %.5dr %10.5dr\n",
         itnum,
         xminp.I[0],
         xminp.I[1],
         fmin);

 } while (status == GSL_CONTINUE && itnum < 100);

gsl_multimin_fdfminimizer_free(s);

/* The End */

record ParamStruct {
  var x0, y0 : real;
  var xscal, yscal : real;
  var fmin : real;
}

proc my_f(v : c_ptr(gsl_vector), params : c_void_ptr) : c_double {
  ref p = params.deref(ParamStruct);

  const x = v.I[0],
        y = v.I[1];
  const retval = p.xscal*(x-p.x0)**2 + p.yscal*(y-p.y0)**2 + p.fmin;
   
  return retval;
}

proc my_df(v : c_ptr(gsl_vector), params : c_void_ptr,
          df : c_ptr(gsl_vector)) {
  ref p = params.deref(ParamStruct);

  const x = v.I[0],
        y = v.I[1];
  df.I[0] = 2*p.xscal*(x-p.x0);
  df.I[1] = 2*p.yscal*(y-p.y0);
}

proc my_fdf(v : c_ptr(gsl_vector), params : c_void_ptr,
            f : c_ptr(c_double), df : c_ptr(gsl_vector)) {
  f.deref() = my_f(v, params);
  my_df(v, params, df);
}
