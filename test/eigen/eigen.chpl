use GSL.Eigen;

var m = new owned GSLMatrix(4,4);
forall (i,j) in m.D do m[i,j] = 1/(i+j+1):real;

var eval = new owned GSLVector(4);
var evec = new owned GSLMatrix(4,4);

var w = gsl_eigen_symmv_alloc(4);

gsl_eigen_symmv(~m, ~eval, ~evec, w);
gsl_eigen_symmv_free(w);

gsl_eigen_symmv_sort(~eval, ~evec, GSL_EIGEN_SORT_ABS_ASC:uint(32));

for i in 0.. #4 {
  writef("eigenvalue = %.8er\n", eval[i]);
  writef("eigenvector = \n");
  for j in 0.. #4 do writef("%.8dr\n",evec[j,i]);
}
