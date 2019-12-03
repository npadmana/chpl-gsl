use GSL.Interpolation;


/* this dataset is taken from
 * J. M. Hyman, Accurate Monotonicity preserving cubic interpolation,
 * SIAM J. Sci. Stat. Comput. 4, 4, 1983. */
var x = [ 7.99, 8.09, 8.19, 8.7, 9.2,
          10.0, 12.0, 15.0, 20.0 ];
var y = [ 0.0, 2.76429e-5, 4.37498e-2,
          0.169183, 0.469428, 0.943740,
          0.998636, 0.999919, 0.999994 ];


var spline_cubic = new Spline(x, y, InterpType.Cubic);
var spline_akima = new Spline(x, y, InterpType.Akima);
var spline_steffen = new Spline(x, y, InterpType.Steffen);

for (xi,yi) in zip(x,y) do writef("%r %r\n",xi,yi);

writef("\n\n");


for i in 0..100 {
  var xi = (1-i/100.0)*x[1] + (i/100.0)*x[x.size];
  var yi_cubic = spline_cubic(xi);
  var yi_akima = spline_akima(xi);
  var yi_steffen = spline_steffen(xi);
  writef("%r %r %r %r\n",xi,yi_cubic, yi_akima, yi_steffen);
}
