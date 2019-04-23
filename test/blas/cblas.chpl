use GSL.Array;

var a = reshape([0.11,0.12,0.13,0.21,0.22,0.23],{0..1,0..2});
var b = reshape([1011.0,1012.0,1021.0,1022.0,1031.0,1032.0],{0..2,0..1});
var c : [0..1,0..1]real;

var lda : c_int = 3,
  ldb : c_int = 2,
  ldc : c_int = 2;

cblas_dgemm(CblasRowMajor : uint(32),
            CblasNoTrans : uint(32),
            CblasNoTrans : uint(32),
            2, 2, 3, 1.0,
            ~a, lda,
            ~b, ldb, 0.0,
            ~c, ldc);

writeln(c);

