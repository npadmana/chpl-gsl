use GSL.SpecFun;

const x = 5.0;
const y = gsl_sf_bessel_J0(x);
writef("J0(%r) = %.18er\n",x,y);
