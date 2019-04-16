use GSL.Interpolation;


/* this dataset is taken from
 * J. M. Hyman, Accurate Monotonicity preserving cubic interpolation,
 * SIAM J. Sci. Stat. Comput. 4, 4, 1983. */
var x = [ 7.99, 8.09, 8.19, 8.7, 9.2,
          10.0, 12.0, 15.0, 20.0 ];
var y = [ 0.0, 2.76429e-5, 4.37498e-2,
          0.169183, 0.469428, 0.943740,
          0.998636, 0.999919, 0.999994 ];
var N = x.size : size_t;

var acc = gsl_interp_accel_alloc();
var spline_cubic = gsl_spline_alloc(gsl_interp_cspline, N);
var spline_akima = gsl_spline_alloc(gsl_interp_akima, N);
var spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N);

gsl_spline_init(spline_cubic, ~x, ~y, N);
gsl_spline_init(spline_akima, ~x, ~y, N);
gsl_spline_init(spline_steffen, ~x, ~y, N);

for (xi,yi) in zip(x,y) do writef("%r %r\n",xi,yi);

writef("\n\n");


for i in 0..100 {
  var xi = (1-i/100.0)*x[1] + (i/100.0)*x[N:int];
  var yi_cubic= gsl_spline_eval(spline_cubic, xi, acc);
  var yi_akima= gsl_spline_eval(spline_akima, xi, acc);
  var yi_steffen= gsl_spline_eval(spline_steffen, xi, acc);
  writef("%r %r %r %r\n",xi,yi_cubic, yi_akima, yi_steffen);
}

gsl_spline_free(spline_cubic);
gsl_spline_free(spline_akima);
gsl_spline_free(spline_steffen);
gsl_interp_accel_free(acc);