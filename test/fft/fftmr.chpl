use GSL.FFT;

const N:size_t=630;

// Domains
const D = {0.. #N};

var data : [D] complex(128);
data = 0.0;
data[0] = 1.0;
forall ii in 1..10 {
  data[ii] = 1.0;
  data[630-ii] = 1.0;
}

for (ii,zz) in zip(D, data) {
  writef("%3i: %13.5er %13.5er\n",ii,zz.re,zz.im);
}
writef("\n");

var wavetable = gsl_fft_complex_wavetable_alloc(N);
var workspace = gsl_fft_complex_workspace_alloc(N);

const nf = wavetable.deref().nf;
const factor = wavetable.deref().factor;
for i in 0.. #nf {
  writef("# factor %i: %u\n", i, factor[i]);
}

gsl_fft_complex_forward((~data):c_ptr(c_double), 1,N,
                        wavetable, workspace);
gsl_fft_complex_wavetable_free(wavetable);
gsl_fft_complex_workspace_free(workspace);

for (ii,zz) in zip(D, data) {
  writef("%3i: %13.5er %13.5er\n",ii,zz.re,zz.im);
}