use GSL.FFT;

const n:size_t=100;
var data : [0.. #n] real;

data = 0.0;
[ii in (n:int/3)..(2*n:int/3-1)] data[ii] = 1.0;

for ii in data.domain do writef("%3i: %13.5er\n", ii, data[ii]);
writef("\n");

var work = gsl_fft_real_workspace_alloc(n);
var rtable = gsl_fft_real_wavetable_alloc(n);

gsl_fft_real_transform(~data, 1, n, rtable, work);
gsl_fft_real_wavetable_free(rtable);

data[11..]=0;

var hc = gsl_fft_halfcomplex_wavetable_alloc(n);

gsl_fft_halfcomplex_inverse(~data, 1, n, hc, work);
gsl_fft_halfcomplex_wavetable_free(hc);

for ii in data.domain do writef("%3i: %13.5er\n", ii, data[ii]);

gsl_fft_real_workspace_free(work);