use GSL.RanDist;

config const x = 2.0;

const P = gsl_cdf_ugaussian_P(x);
const Q = gsl_cdf_ugaussian_Q(x);
const x1 = gsl_cdf_ugaussian_Pinv(P);
const x2 = gsl_cdf_ugaussian_Qinv(Q);

writef("prob(x < %8.6dr) = %8.6dr\n",x,P);
writef("prob(x > %8.6dr) = %8.6dr\n",x,Q);
writef("Pinv(%8.6dr) = %8.6dr\n",P,x1);
writef("Qinv(%8.6dr) = %8.6dr\n",Q,x2);