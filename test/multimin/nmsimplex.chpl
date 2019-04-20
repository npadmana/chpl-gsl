/* Nelder Mead minimization example */

use GSL.Minimize;

// Starting point
var x : [0.. #2] real = [5.0, 7.0];

// Step sizes
var ss : [0.. #2] real = [1.0, 1.0];

// Parameters of parabaloid
var pp = new ParamStruct(1.0, 2.0, 10.0, 20.0, 30.0);

// Set up user function and minimizer
var minex_func : gsl_multimin_function;
minex_func.n = 2;
minex_func.f = c_ptrTo(my_f);
minex_func.params = c_ptrTo(pp);

var s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, 2);
// Get GSL vector views on arrays
var xview = GSLView(x);
var ssview = GSLView(ss);
gsl_multimin_fminimizer_set(s, ~minex_func,
                            ~xview,
                            ~ssview);

// Iterations
var status : c_int;
var itnum = 0;

do {
  itnum += 1;
  status = gsl_multimin_fminimizer_iterate(s);

  if status then break;

  var size = gsl_multimin_fminimizer_size(s);
  status = gsl_multimin_test_size(size, 1.0e-2);

  if status==GSL_SUCCESS then writef("converged to minimum at\n");
  var xminp = gsl_vector_ptr(gsl_multimin_fminimizer_x(s), 0);
  var fmin = gsl_multimin_fminimizer_minimum(s);
  var fsize = gsl_multimin_fminimizer_size(s);
  writef("%5i %10.3er %10.3er f() = %7.3dr size = %.3dr\n",
         itnum,
         xminp[0],
         xminp[1],
         fmin,
         fsize);

 } while (status == GSL_CONTINUE && itnum < 100);

gsl_multimin_fminimizer_free(s);

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
