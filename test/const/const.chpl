use GSL.Constants;

const c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
const au = GSL_CONST_MKSA_ASTRONOMICAL_UNIT;
const minutes = GSL_CONST_MKSA_MINUTE;

const r_earth = 1.0*au;
const r_mars = 1.52*au;

writef("light travel time from Earth to Mars:\n");
writef("minimum = %.1dr minutes\n", (r_mars-r_earth)/(c*minutes));
writef("maximum = %.1dr minutes\n", (r_mars+r_earth)/(c*minutes));