use GSL.IEEE;

extern proc setenv(name: c_string, val : c_string, ovride : c_int) : c_int;

config const choice=1;
select choice {
  when 1 do setenv("GSL_IEEE_MODE".c_str(),
                   "round-to-nearest".c_str(),
                   1);
  when 2 do setenv("GSL_IEEE_MODE".c_str(),
                   "round-down",
                   1);
}

gsl_ieee_env_setup();
dosum();

proc dosum() {

  var x : real(64) = 1.0,
    sum : real(64) = 0.0,
    oldsum : real(64) = 0.0;

  var i = 0;

  do {
    i += 1;
    oldsum = sum;
    sum += x;
    x = x/i;

    writef("i=%2i sum=%.18dr error=%r\n",
           i, sum, sum - M_E);
    if i > 30 then break;
  } while (sum != oldsum);
}
