use GSL.Statistics;

var data : [0..#5]real = [17.2, 18.1, 16.5, 18.3, 12.6];

var mean     = gsl_stats_mean(~data, 1, 5);
var variance = gsl_stats_variance(~data, 1, 5);
var largest  = gsl_stats_max(~data, 1, 5);
var smallest = gsl_stats_min(~data, 1, 5);

writef ("The dataset is %r, %r, %r, %r, %r\n",
        data[0], data[1], data[2], data[3], data[4]);

writef ("The sample mean is %r\n", mean);
writef ("The estimated variance is %r\n", variance);
writef ("The largest value is %r\n", largest);
writef ("The smallest value is %r\n", smallest);
