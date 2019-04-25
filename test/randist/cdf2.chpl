use GSL.RanDist;

config const x = 2.0;
var gauss = new Gaussian(1.0);

const P = gauss.P(x);
const Q = gauss.Q(x);
const x1 = gauss.Pinv(P);
const x2 = gauss.Qinv(Q);

writef("prob(x < %8.6dr) = %8.6dr\n",x,P);
writef("prob(x > %8.6dr) = %8.6dr\n",x,Q);
writef("Pinv(%8.6dr) = %8.6dr\n",P,x1);
writef("Qinv(%8.6dr) = %8.6dr\n",Q,x2);