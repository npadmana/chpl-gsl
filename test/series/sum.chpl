use GSL.SeriesAccel;

config const N : size_t = 20;

const zeta_2 = pi**2/6.0;
var t : [0.. #20] real;
var sum : real = 0.0;
var sum_accel, err : real;

for n in t.domain {
  var np1 = n + 1.0;
  t[n] = 1.0/np1**2;
  sum += t[n];
}
var w = gsl_sum_levin_u_alloc(N);
gsl_sum_levin_u_accel(~t, N, w, ~sum_accel,~err);

writef ("term-by-term sum = % .16dr using %i terms\n", 
        sum, N);
writef ("term-by-term sum = % .16dr using %i terms\n", 
        (w.deref()).sum_plain, (w.deref()).terms_used);

writef ("exact value      = % .16dr\n", zeta_2);
writef ("accelerated sum  = % .16dr using %i terms\n", 
        sum_accel, (w.deref()).terms_used);

writef ("estimated error  = % .16dr\n", err);
writef ("actual error     = % .16dr\n", 
        sum_accel - zeta_2);
gsl_sum_levin_u_free(w);