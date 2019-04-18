use GSL.Eigen;

// Set up the matrix
var m = new owned GSLMatrix(4,4);
{
  var x:[0..3]real = [-1.0,-2.0,3.0,4.0];
  forall (i,j) in m.D {
    m(i,j) = x[i]**(3-j);
  }
}

// Set up matrices and vectors for the eigen problem
var eval = new owned GSLVector(complex(128),4);
var evec = new owned GSLMatrix(complex(128),4,4);

{
  var w = gsl_eigen_nonsymmv_alloc(4);
  gsl_eigen_nonsymmv(~m, ~eval, ~evec, w);
  gsl_eigen_nonsymmv_free(w);
  gsl_eigen_nonsymmv_sort(~eval, ~evec, GSL_EIGEN_SORT_ABS_DESC:uint(32));
}


for i in 0.. #4 {
  writef("eigenvalue = %.8dr + %.8dri\n", eval[i].re, eval[i].im);
  writef("eigenvector = \n");
  for j in 0.. #4 do writef("%.8dr + %.8dri\n",evec[j,i].re, evec[j,i].im);
}