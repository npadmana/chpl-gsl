use GSL.SpecFun;

const x = 5.0;
const expected = -0.17759677131433830434739701;
  
const y = gsl_sf_bessel_J0 (x);

writef ("J0(5.0) = %.18dr\n", y);
writef ("exact   = %.18dr\n", expected);