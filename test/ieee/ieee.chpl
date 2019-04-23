use GSL.IEEE;

var f : real(32) = (1.0/3.0):real(32);
var d : real(64) = 1.0/3.0;
var fd = f : real(64);

writef(" f=");
gsl_ieee_printf_float(c_ptrTo(f));
writef("\n");

writef("fd=");
gsl_ieee_printf_double(c_ptrTo(fd));
writef("\n");

writef(" d=");
gsl_ieee_printf_double(c_ptrTo(d));
writef("\n");