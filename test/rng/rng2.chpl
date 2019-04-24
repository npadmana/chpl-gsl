use GSL.RNG;

config const RNGType=RNGAlgorithms.MT19937;
config const iseed : uint(64) = 0;

var r = new owned Random(RNGType, iseed);

writef("generator type: %s\n", r.name);
writef("seed = %i\n", r.seed);
writef("first value = %i\n", r.getRaw());



