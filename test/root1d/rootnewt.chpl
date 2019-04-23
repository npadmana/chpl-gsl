use GSL.Roots;
use demo;

const r_expected = sqrt(5.0);
var params = new quadratic_params(1.0, 0.0, -5.0);

var FDF = mkGSLFunc(c_ptrTo(quadratic),
                    c_ptrTo(quadratic_deriv),
                    c_ptrTo(quadratic_fdf),
                    params);

var s = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
var x = 5.0;
gsl_root_fdfsolver_set(s, ~FDF, x);

writef("using %s method\n",gsl_root_fdfsolver_name(s):string);

writef("%-5s %10s %10s %10s\n",
       "iter", "root", "err", "err(est)");

var status : c_int = GSL_CONTINUE;
var itnum = 0, max_iter=100;
do {
  itnum +=1;
  status = gsl_root_fdfsolver_iterate(s);
  var x0 = x;
  x = gsl_root_fdfsolver_root(s);
  status = gsl_root_test_delta(x, x0, 0, 1.0e-3);
  if status==GSL_SUCCESS then writef("Converged:\n");
  writef("%5i %10.7dr %+10.7dr %10.7dr\n",
         itnum, x, x - r_expected, x - x0);
 }
while (status == GSL_CONTINUE && itnum < max_iter);

gsl_root_fdfsolver_free(s);


