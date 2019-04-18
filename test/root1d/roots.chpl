use GSL.Root1D;
use demo;

const r_expected = sqrt(5.0);
var params = new quadratic_params(1.0, 0.0, -5.0);

var F = mkGSLFunc(c_ptrTo(quadratic), params);

var s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
var x_lo : real = 0.0,
  x_hi : real = 5.0;
gsl_root_fsolver_set(s, ~F, x_lo, x_hi);

writef("using %s method\n",gsl_root_fsolver_name(s):string);

writef("%5s [%9s, %9s] %9s %10s %9s\n",
       "iter", "lower", "upper", "root", 
       "err", "err(est)");

var status : c_int = GSL_CONTINUE;
var itnum = 0, max_iter=100;
do {
  itnum +=1;
  status = gsl_root_fsolver_iterate(s);
  var r = gsl_root_fsolver_root(s);
  x_lo = gsl_root_fsolver_x_lower(s);
  x_hi = gsl_root_fsolver_x_upper(s);
  status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);
  if status==GSL_SUCCESS then writef("Converged:\n");

  writef("%5i [%.7dr, %.7dr] %.7dr %+.7dr %.7dr\n",
         itnum, x_lo, x_hi,
         r, r - r_expected, 
         x_hi - x_lo);
 }
while (status == GSL_CONTINUE && itnum < max_iter);

gsl_root_fsolver_free(s);


