/* A simple linear fitting example,
   including weights.

*/
use GSL.LinearFit;

var n : uint(64) = 4;

var x = [1970.0, 1980.0, 1990.0, 2000.0];
var y = [12.0, 11.0, 14.0, 13.0];
var w = [0.1, 0.2, 0.3, 0.4];

// Define variables for fit
var c0, c1, cov00, cov01, cov11, chisq : real;

gsl_fit_wlinear(~x, 1, ~w, 1, ~y, 1, n,
                ~c0, ~c1, ~cov00,
                ~cov01, ~cov11,
                ~chisq);

writef("# best fit: Y = %r + %r X\n", c0, c1);
writef ("# covariance matrix:\n");
writef ("# [ %r, %r\n#   %r, %r]\n", 
        cov00, cov01, cov01, cov11);
writef ("# chisq = %r\n", chisq);

// Write out the data
for (xi,yi,ei) in zip(x,y,w) do writef("data: %r %r %r\n",xi,yi,1.0/sqrt(ei));
writef("\n");

for i in -30..129 {
  var xf = x[1] + (x[n:int]-x[1]) * (i/100.0);
  var yf, yf_err : real;

  gsl_fit_linear_est(xf, c0, c1, cov00, cov01, cov11,
                     ~yf, ~yf_err);

  writef ("fit: %r %r\n", xf, yf);
  writef ("hi : %r %r\n", xf, yf + yf_err);
  writef ("lo : %r %r\n", xf, yf - yf_err);
}