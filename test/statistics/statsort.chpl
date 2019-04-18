use GSL.Statistics;
use GSL.Sorting;

var data : [0..#5]real = [17.2, 18.1, 16.5, 18.3, 12.6];

writef ("Original dataset:  %r, %r, %r, %r, %r\n",
        data[0], data[1], data[2], data[3], data[4]);

gsl_sort(~data, 1, 5);
writef ("Sorted dataset: %r, %r, %r, %r, %r\n",
        data[0], data[1], data[2], data[3], data[4]);

const median = gsl_stats_median_from_sorted_data(~data,1,5);
const upperq = gsl_stats_quantile_from_sorted_data(~data,1,5,0.75);
const lowerq = gsl_stats_quantile_from_sorted_data(~data,1,5,0.25);

writef ("The median is %r\n", median);
writef ("The upper quartile is %r\n", upperq);
writef ("The lower quartile is %r\n", lowerq);
