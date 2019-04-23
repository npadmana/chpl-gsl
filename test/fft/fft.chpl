use GSL.FFT;

const N=128;

// Domains
const D = {0.. #N};

var data : [D] complex(128);
data = 0.0;
data[0] = 1.0;
forall ii in 1..10 {
  data[ii] = 1.0;
  data[128-ii] = 1.0;
}

for (ii,zz) in zip(D, data) {
  writef("%3i %13.5er %13.5er\n",ii,zz.re,zz.im);
}
writef("\n\n");

gsl_fft_complex_radix2_forward((~data):c_ptr(c_double),1,128);
data /= sqrt(128.0);

for (ii,zz) in zip(D, data) {
  writef("%3i %13.5er %13.5er\n",ii,zz.re,zz.im);
}