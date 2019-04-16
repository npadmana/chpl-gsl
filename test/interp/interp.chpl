use GSL.Interpolation;

var x,y :[0.. #10]real;

writef("#m=0,S=17\n");
forall ii in x.domain {
  x[ii] = ii+0.5*sin(ii);
  y[ii] = ii + cos(ii*ii);
}
for (xi,yi) in zip(x,y) do writef("%r %r\n",xi,yi);

writef("#m=1,S=0\n");

var acc = gsl_interp_accel_alloc();
var spline = gsl_spline_alloc(gsl_interp_cspline, 10);
gsl_spline_init(spline, ~x, ~y, 10);

{
  var xi = x[0];
  while (xi < x[9]) {
    var yi = gsl_spline_eval(spline, xi, acc);
    writef("%r %r\n",xi,yi);
    xi += 0.01;
  }
}

gsl_spline_free(spline);
gsl_interp_accel_free(acc);

