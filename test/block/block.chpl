// Use GSL blocks.
// This is a direct translation of the code
// but it is quite unlikely we will need to use
// this code much.

// We also don't print out the block address, since
// that changes and is not a good test.

use GSL.Array;

var b = gsl_block_alloc(100);
writef("length of block = %u\n", (b.deref()).size);
gsl_block_free(b);
