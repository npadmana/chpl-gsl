WARNING_FLAGS=--ccflags="-Wno-incompatible-pointer-types"


nmsimplex.x : nmsimplex.chpl
	chpl -o nmsimplex.x nmsimplex.chpl ${GSL_FLAGS} ${WARNING_FLAGS}


cleanest:
	-rm *.x 
