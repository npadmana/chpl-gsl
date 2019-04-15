use GSL.SpecFun;

const x = 5.0;
const expected = -0.17759677131433830434739701;

var result : gsl_sf_result;

var status = gsl_sf_bessel_J0_e(x, ~result);

writef("status  = %s\n", gsl_strerror(status):string);
writef("J0(5.0) = %.18dr\n", result.val);
writef("      +/-  %.18dr\n", result.err);
writef("exact   = %.18dr\n", expected);