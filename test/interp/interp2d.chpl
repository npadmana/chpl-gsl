use GSL.Interpolation;

const N : size_t = 100;

var xa = [0.0, 1.0];
var ya = [0.0, 1.0];
const nx = xa.size: size_t;
const ny = ya.size: size_t;
var za : [0..#nx, 0..#ny] real; // Use a 2D array

var spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
var xacc = gsl_interp_accel_alloc();
var yacc = gsl_interp_accel_alloc();

za[0,0]=0.0;
za[0,1]=1.0;
za[1,1]=0.5;
za[1,0]=1.0;

gsl_spline2d_init(spline, ~xa, ~ya, ~za, nx,ny);

for i in 0.. #N {
  var xi = i/(N-1.0);
  for j in 0.. #N {
    var yj = j/(N-1.0);
    var zij = gsl_spline2d_eval(spline, xi, yj, xacc, yacc);

    writef("%dr %dr %dr\n",xi,yj,zij);
  }
  writef("\n");
}

gsl_spline2d_free(spline);
gsl_interp_accel_free(xacc);
gsl_interp_accel_free(yacc);