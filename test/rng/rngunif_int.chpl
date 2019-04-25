use GSL.RNG;

config const RNGType=RNGAlgorithms.MT19937;
config const iseed : uint(64) = 0;

var r = new Random(RNGType, iseed);
for i in 1..10 do writef("%i\n", r.getInt(10));



