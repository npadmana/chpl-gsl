use GSL.Interpolation;

var x = [0.00, 0.10,  0.27, 0.30];
var y = [0.15, 0.70, -0.10, 0.15];
// Note that y[0]==y[3] for periodic data

writef("#m=0,S=5\n");
for (xi,yi) in zip(x,y) do writef("%r %r\n",xi,yi);


var acc = gsl_interp_accel_alloc();
var spline = gsl_spline_alloc(gsl_interp_cspline_periodic, x.size : size_t);
gsl_spline_init(spline, ~x, ~y, x.size : size_t);

writef("#m=1,S=0\n");
for i in 0..100 {
  var xi = (1-i/100.0)*x[1] + (i/100.0)*x[4];
  var yi = gsl_spline_eval(spline, xi, acc);
  writef("%r %r\n",xi,yi);
}

gsl_spline_free(spline);
gsl_interp_accel_free(acc);