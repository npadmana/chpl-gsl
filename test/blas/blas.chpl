use GSL.Array;

var a = reshape([0.11,0.12,0.13,0.21,0.22,0.23],{0..1,0..2});
var b = reshape([1011.0,1012.0,1021.0,1022.0,1031.0,1032.0],{0..2,0..1});
var c : [0..1,0..1]real;

var A = GSLView(a);
var B = GSLView(b);
var C = GSLView(c);

var ret = gsl_blas_dgemm(CblasNoTrans:uint(32),
                         CblasNoTrans:uint(32),
                         1.0,
                         ~A, ~B, 0.0, ~C);

writeln(c);
