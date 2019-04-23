/* Read in a signal, perform a DWT,
   keep only the 20 largest coefficients,
   and then inverse DWT. */

use GSL.DWT;
use GSL.Sorting;

const N:uint(64)=256; // File size
const NC:uint(64)=20; // Number of non-zero components.
const fn="ecg.dat"; // Name of the file

var orig_data, data: [0.. #N]real;
// Read in data
var ff = openreader(fn);
for x in orig_data do ff.read(x);
ff.close();
data = orig_data;

var w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
var work = gsl_wavelet_workspace_alloc(N);

gsl_wavelet_transform_forward(w, ~data, 1, N, work);

// Get the absolute coefficients 
var abscoeff = abs(data);
var p : [0.. #N] uint(64);
gsl_sort_index(~p, ~abscoeff, 1, N);
[ii in 0.. #(N-NC)] data[p[ii]:int]=0.0;

// Inverse wavelet transform
gsl_wavelet_transform_inverse(w, ~data, 1, N, work);

for (x,y) in zip(orig_data, data) do writef("%.8dr %.8dr\n",x,y);

gsl_wavelet_free(w);
gsl_wavelet_workspace_free(work);





