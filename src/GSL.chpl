/*
  This is the top-level GSL module. Most of the GSL functionality
  is included in a series of submodules, which are not automatically
  included.

  Instead of explicit declarations, this makes heavy use of extern blocks. The
  Chapel code here is restricted to convenience functions and workarounds.

  The mapping from sub-modules to GSL chapters is :

  * :mod:`Array` : Vectors and Matrices, BLAS, CBLAS
  * :mod:`BSpline` : Basis Splines
  * :mod:`Chebyshev` : Chebyshev Approximations
  * :mod:`Combinations` : Combinations
  * :mod:`Complex` : Complex numbers
  * :mod:`Constants` : Physical Constants
  * :mod:`Deriv` : Numerical Differentiation
  * :mod:`DHT` : Discrete Hankel Transforms
  * :mod:`DWT` : Discrete Wavelet Transforms
  * :mod:`Eigen` : Eigensystems
  * :mod:`Fit` : Linear and Nonlinear Least-Squares Fitting
  * :mod:`FFT` : Fast Fourier Transforms
  * :mod:`Filter` : Digital Filtering
  * :mod:`Histograms` : Histograms
  * :mod:`IEEE` : IEEE floating-point arithmetic
  * :mod:`Integration` : Numerical Integration
  * :mod:`Interpolation` : Interpolation
  * :mod:`LinearAlgebra` : Linear Algebra
  * :mod:`Minimize` : One- and Multi-dimensional Minimization
  * :mod:`MonteCarlo` : Monte Carlo Integration
  * :mod:`Multiset` : Multisets
  * :mod:`Ntuple` : N-tuples 
  * :mod:`ODE` : Ordinary Differential Equations
  * :mod:`Permutations` : Permutations
  * :mod:`Polynomials` : Polynomials
  * :mod:`QRNG` : Quasi-random sequences
  * :mod:`RNG` : Random number generation
  * :mod:`RanDist` : Random number distributions
  * :mod:`Roots` : One- and Multi-dimensional root finding
  * :mod:`SeriesAccel` : Series Acceleration
  * :mod:`SimAnneal` : Simulated Annealing
  * :mod:`Sorting` : Sorting
  * :mod:`Sparse` : Sparse Matrices, BLAS and Linear Algebra
  * :mod:`SpecFun` : Special Functions
  * :mod:`Statistics` : Statistics, Running Statistics, Moving Window Statistics
  
*/
module GSL {

  use Common;

  module Common {

    extern record gsl_function {
      var function : c_fn_ptr;
      var params : c_void_ptr;
    }

    extern record gsl_function_fdf {
      var f : c_fn_ptr;
      var df : c_fn_ptr;
      var fdf : c_fn_ptr;
      var params : c_void_ptr;
    }

    extern {
      #include "gsl/gsl_errno.h"
      #include "gsl/gsl_math.h"
    }


    /*
      Dereference a C void pointer to an arbitrary type.

      This is useful when converting parameters in callbacks.
    */
    inline proc c_void_ptr.deref(type t) ref {
      return (this:c_ptr(t)).deref();
    }

    /*
      Overload the ``~`` operator to return a C pointer to an array.
    */
    inline proc ~(ref x : []?T): c_ptr(T) {
      return c_ptrTo(x);
    }

    /*
      Overload the ``~`` operator to return a C pointer.

      .. warning::

      This will not work for ``bool`` types.
    */
    inline proc ~(ref x : ?T) : c_ptr(T)
      where !isArrayType(T)
    {
      return c_ptrTo(x);
    }

    /* Take a function pointer and parameters and wrap it
       into a ``gsl_function`` construct.
    */
    proc mkGSLFunc(fptr, ref params) {
      var F : gsl_function;
      F.function = fptr : c_fn_ptr;
      F.params = c_ptrTo(params);
      return F;
    }

    /* Take a function pointer and wrap it into a ``gsl_function`` construct.

       The ``params`` part is set  to void in this case.
    */
    proc mkGSLFunc(fptr) {
      var F : gsl_function;
      F.function = fptr : c_fn_ptr;
      F.params = nil : c_void_ptr;
      return F;
    }

    /* Take a function pointer and parameters and wrap it
       into a ``gsl_function_fdf`` construct.
    */
    proc mkGSLFunc(f,df,fdf, ref params) {
      var F : gsl_function_fdf;
      F.f = f : c_fn_ptr;
      F.df = df : c_fn_ptr;
      F.fdf = fdf : c_fn_ptr;
      F.params = c_ptrTo(params);
      return F;
    }
    
    /* Take a function pointer and parameters and wrap it
       into a ``gsl_function_fdf`` construct.

       The ``params`` part is set  to void in this case.
    */
    proc mkGSLFunc(f,df,fdf) {
      var F : gsl_function_fdf;
      F.f = f : c_fn_ptr;
      F.df = df : c_fn_ptr;
      F.fdf = fdf : c_fn_ptr;
      F.params = nil : c_void_ptr;
      return F;
    }

  }

  /* GSL Matrices and Vectors

     This includes support for the ``gsl_matrix`` and ``gsl_vector``
     structures. Note that, in the process of including the headers
     for these, we end up also including headers for blocks and
     other structures, so we don't explicitly include these.

     The code also includes helper routines/classes for these structures
     to make them easier to deal with in Chapel.

     I expect this modules to be automatically included in a number of
     different places below.

     .. note::

        There are a number of improvements still possible.
        * Iterators for the vector and matrix classes
        * The matrix classes should support more than just ``real(64)``

  */
  module Array {

    // TODO : Extern blocks don't translate these well.
    extern record gsl_matrix {
      var data : c_ptr(c_double);
    }
    extern record gsl_vector {
      var data : c_ptr(c_double);
    }
    extern record gsl_matrix_complex {
      var data : c_ptr(c_double);
    }
    extern record gsl_vector_complex {
      var data : c_ptr(c_double);
    }
    extern record gsl_matrix_int {
      var data : c_ptr(c_int);
    }
    extern record gsl_vector_int {
      var data : c_ptr(c_int);
    }
    extern record gsl_vector_view {
      var vector : gsl_vector;
    }
    extern record gsl_matrix_view {
      var matrix : gsl_matrix;
    }

    extern {
      #include "gsl/gsl_vector.h"
      #include "gsl/gsl_matrix.h"
      #include "gsl/gsl_blas.h"
      #include "gsl/gsl_cblas.h"

      // Unfortunately, GSL returns a slightly different type and so
      // we write some helper routines here to clean things up.
      static gsl_matrix_view np_gsl_matrix_view_array(double* arr, size_t s1, size_t s2) {
        return gsl_matrix_view_array(arr,s1,s2);
      }

      static gsl_vector_view np_gsl_vector_view_array(double* arr, size_t s1) {
        return gsl_vector_view_array(arr,s1);
      }
    }


    /* Create a GSL view on a 2D array */
    proc GSLView(arr : [?D] real) : gsl_matrix_view
      where (D.rank==2) {
      return np_gsl_matrix_view_array(c_ptrTo(arr), D.dim(1).size:size_t,
                                   D.dim(2).size:size_t);
    }

    /* Create a GSL view on a 1D array */
    proc GSLView(arr : [?D] real) 
      where (D.rank==1) {
      return np_gsl_vector_view_array(c_ptrTo(arr), D.dim(1).size:size_t);
    }

    /* Return a pointer to a GSL matrix view */
    inline proc ~(ref v : gsl_matrix_view) {
      return c_ptrTo(v.matrix);
    }

    /* Return a pointer to a GSL vector view */
    inline proc ~(ref v : gsl_vector_view) {
      return c_ptrTo(v.vector);
    }

    /* Wrap the GSL Matrix type */
    class GSLMatrix {
      type eltType;
      type dataType;
      type gslMatType;
      var p : gslMatType; // Pointer to GSL Matrix
      var pdata : c_ptr(dataType); // Pointer to the actual data
      var D : domain(2);

      proc init(m,n) where (isIntegral(m) && isIntegral(n)) {
        eltType = c_double;
        dataType = c_double;
        gslMatType = c_ptr(gsl_matrix);
        this.complete();
        p = gsl_matrix_alloc(m:size_t, n:size_t);
        pdata = (p.deref()).data;
        D = {0.. #m, 0.. #n};
      }

      proc init(type t, m, n) where (isIntegral(m) && isIntegral(n)) {
        eltType = t;
        select t {
            when real(64) {
              dataType = c_double;
              gslMatType = c_ptr(gsl_matrix);
            }
            when c_double {
              dataType = c_double;
              gslMatType = c_ptr(gsl_matrix);
            }
            when c_int {
              dataType = c_int;
              gslMatType = c_ptr(gsl_matrix_int);
            }
            when complex(128) {
              dataType = c_double; // Complex data is stored as a list of doubles
              gslMatType = c_ptr(gsl_matrix_complex);
            }
            otherwise {
              dataType = bool;
              gslMatType = bool;
              compilerError("Unimplemented type for GSLMatrix : ",t:string);
            }
          }
        this.complete();
        select t {
            when real(64) do p = gsl_matrix_alloc(m:size_t, n:size_t);
            when c_double do p = gsl_matrix_alloc(m:size_t, n:size_t);
            when c_int do p = gsl_matrix_int_alloc(m:size_t, n:size_t);
            when complex(128) do p = gsl_matrix_complex_alloc(m:size_t, n:size_t);
          }
        pdata = (p.deref()).data;
        D = {0.. #m, 0.. #n};
      }

      proc deinit() {
        select eltType {
            when real(64) do gsl_matrix_free(p);
            when c_double do gsl_matrix_free(p);
            when c_int do gsl_matrix_int_free(p);
            when complex(128) do gsl_matrix_complex_free(p);
          }
      }

      proc access(i : size_t, j : size_t) ref {
        select eltType {
            when real(64) do return (gsl_matrix_ptr(p,i,j)).deref();
            when c_double do return (gsl_matrix_ptr(p,i,j)).deref();
            when c_int do return (gsl_matrix_int_ptr(p,i,j)).deref();
            when complex(128) do return (gsl_matrix_complex_ptr(p,i,j):c_ptr(complex(128))).deref();
          }
      }

      proc this(i,j) ref where (isIntegral(i) && isIntegral(j)) {
        return access(i:size_t, j:size_t);
      }
    }

    /* Return a pointer to a GSLMatrix. In this case, we return the pointer to
       ``gsl_matrix`` which is what one normally needs.
    */
    inline proc ~(m : GSLMatrix) {
      return m.p;
    }

    /* Wrap the GSL Vector type */
    class GSLVector {
      type eltType;
      type dataType;
      type gslVecType;
      var p : gslVecType; // Pointer to GSL Vector
      var pdata : c_ptr(dataType); // Pointer to the actual data
      var D : domain(1);

      proc init(n) where isIntegral(n) {
        eltType = c_double;
        dataType = c_double;
        gslVecType = c_ptr(gsl_vector);
        this.complete();
        p = gsl_vector_alloc(n:size_t);
        pdata = (p.deref()).data;
        D = {0.. #n};
      }

      proc init(type t, n) where isIntegral(n) {
        eltType = t;
        select t {
            when real(64) {
              dataType = c_double;
              gslVecType = c_ptr(gsl_vector);
            }
            when c_double {
              dataType = c_double;
              gslVecType = c_ptr(gsl_vector);
            }
            when c_int {
              dataType = c_int;
              gslVecType = c_ptr(gsl_vector_int);
            }
            when complex(128) {
              dataType = c_double; // Complex data is stored as a list of doubles
              gslVecType = c_ptr(gsl_vector_complex);
            }
            otherwise {
              dataType = bool;
              gslVecType = bool;
              compilerError("Unimplemented type for GSLVector : ",t:string);
            }
          }
        this.complete();
        select t {
            when real(64) do p = gsl_vector_alloc(n:size_t);
            when c_double do p = gsl_vector_alloc(n:size_t);
            when c_int do p = gsl_vector_int_alloc(n:size_t);
            when complex(128) do p = gsl_vector_complex_alloc(n:size_t);
          }
        pdata = (p.deref()).data;
        D = {0.. #n};
      }

      proc deinit() {
        select eltType {
            when real(64) do gsl_vector_free(p);
            when c_double do gsl_vector_free(p);
            when c_int do gsl_vector_int_free(p);
            when complex(128) do gsl_vector_complex_free(p);
          }
      }

      proc access(i : size_t) ref {
        select eltType {
            when real(64) do return (gsl_vector_ptr(p,i)).deref();
            when c_double do return (gsl_vector_ptr(p,i)).deref();
            when c_int do return (gsl_vector_int_ptr(p,i)).deref();
            when complex(128) do return (gsl_vector_complex_ptr(p,i):c_ptr(complex(128))).deref();
          }
      }

      proc this(i) ref where isIntegral(i) {
        return access(i:size_t);
      }
    }

    /* Return a pointer to a GSLVector. In this case, we return the pointer to
       ``gsl_vector`` which is what one normally needs.

       Note that, since ``v`` is a class, we don't need to decorate it with a ``ref``???.
       The issue appears to be the fact that GSLVector is generic.
    */
    inline proc ~(v : GSLVector) {
      return v.p;
    }

    /* Simple accessor functions for GSL vectors

       .. note::

          Why does ``this`` not work here?
       
     */
    proc (c_ptr(gsl_vector)).I(i) ref {
      return (gsl_vector_ptr(this,i:size_t)).deref();
    }

    /* Simple accessor functions for GSL vectors */
    proc (c_ptr(gsl_matrix)).I(i,j) ref {
      return (gsl_matrix_ptr(this,i:size_t, j:size_t)).deref();
    }

  }

  /* Basis Splines

     Includes ``gsl_bspline.h``

     This automatically includes ``Array`` and ``Common`` as well.
  */
  module BSpline {
    use Common;
    use Array;

    extern {
      #include "gsl/gsl_bspline.h"
    }
  }

  /* Combinations

     Includes ``gsl_combinations.h``
  */
  module Combinations {
    extern {
      #include "gsl/gsl_combination.h"
    }
  }

  /* Chebyshev approximations

     Includes ``gsl_chebyshev.h``
  */
  module Chebyshev {
    extern {
      #include "gsl/gsl_chebyshev.h"
    }
  }


  /* Support for GSL complex numbers.

     This includes both the ``gsl_complex.h`` and ``gsl_complex_math.h``
     libraries.

     .. note::

     We explicitly define the ``gsl_complex`` and ``gsl_complex_float``
     types, to avoid using the macros that GSL uses. These decisions
     may get revisited later.

  */
  module Complex {

    extern record gsl_complex {
      var dat : 2*c_double;
    }
    //extern type gsl_complex = complex(128);

    extern record gsl_complex_float {
      var dat : 2*c_float;
    }

    extern {
      #include "gsl/gsl_complex.h"
      #include "gsl/gsl_complex_math.h"
    }

    /* Convert from Chapel complex to gsl_complex */
    inline proc Z_(x:complex(128)):gsl_complex {
      return new gsl_complex((x.re,x.im));
    }

    /* Convert from gsl_complex to Chapel complex */
    inline proc Z_(x:gsl_complex):complex(128) {
      return x.dat:complex(128);
    }

    /* Convert from Chapel complex to gsl_complex */
    inline proc Z_(x:complex(64)):gsl_complex_float {
      return new gsl_complex_float((x.re,x.im));
    }

    /* Convert from gsl_complex to Chapel complex */
    inline proc Z_(x:gsl_complex_float):complex(64) {
      return x.dat:complex(64);
    }

    // End module complex.
  }


  /* Physical Constants

     Includes ``gsl_const_{mksa,cgsm,num}.h``
  */
  module Constants {
    extern {
      #include "gsl/gsl_const_mksa.h"
      #include "gsl/gsl_const_cgsm.h"
      #include "gsl/gsl_const_num.h"
    }
  }

  /* Numerical Differentiation.

     Based on ``gsl_deriv.h``.
  */
  module Deriv {
    use Common;

    extern {
      #include "gsl/gsl_deriv.h"
    }
  }

  /* Hankel Transforms

     Includes ``gsl_dht.h``
  */
  module DHT {
    extern {
      #include "gsl/gsl_dht.h"
    }
  }


  /* Wavelet Transforms

     Includes ``gsl_wavelet.h`` and ``gsl_wavelet2d.h``
  */
  module DWT {
    extern {
      #include "gsl/gsl_wavelet.h"
      #include "gsl/gsl_wavelet2d.h"
    }
  }


  /* Eigensystems

     Based on ``gsl_eigen.h``
  */
  module Eigen {
    use Common;
    use Array;

    extern {
      #include "gsl/gsl_eigen.h"
    }

  }

  /* Fast Fourier Transforms

     Includes :
     * ``gsl_fft_complex.h``
     * ``gsl_fft_real.h``
     * ``gsl_fft_halfcomplex.h``

  */
  module FFT {
    use Complex;

    extern {
      #include "gsl/gsl_fft_complex.h"
      #include "gsl/gsl_fft_real.h"
      #include "gsl/gsl_fft_halfcomplex.h"
    }
  }


  /* Digital Filtering

     Includes ``gsl_filter.h``
  */
  module Filter {
    use Array;

    extern {
      #include "gsl/gsl_filter.h"
    }
  }

  /* Histograms

     Includes ``gsl_histogram.h``, ``gsl_histogram2d.h``
  */
  module Histograms {
    extern {
      #include "gsl/gsl_histogram.h"
      #include "gsl/gsl_histogram2d.h"
    }
  }

  /* IEEE Arithmetic

   */
  module IEEE {
    extern {
      #include "gsl/gsl_ieee_utils.h"
    }
  }



  /* Numerical integration

     Based on ``gsl_integration.h`` 
  */
  module Integration {
    use Common;

    extern {
      #include "gsl/gsl_integration.h"
    }
  }


  /* Interpolation

     Based on ``gsl_interp.h``, ``gsl_spline.h`` and the
     2D versions of the same.

  */
  module Interpolation {
    extern {
      #include "gsl/gsl_interp.h"
      #include "gsl/gsl_spline.h"
      #include "gsl/gsl_interp2d.h"
      #include "gsl/gsl_spline2d.h"
    }
  }


  /* Linear Algebra

     Based on
     * ``gsl_linalg.h``

  */
  module LinearAlgebra {
    use Common;
    use Array;

    extern {
      #include "gsl/gsl_linalg.h"
    }
  }

  /* Least-Squares Fitting

     Based on ``gsl_fit.h``, ``gsl_multifit.h`` and
     ``gsl_multilarge.h``, ``gsl_multifit_nlinear.h``
     and ``gsl_multilarge_nlinear.h``
   */
  module Fit {
    use Common;
    use Array;

    extern {
      #include "gsl/gsl_fit.h"
      #include "gsl/gsl_multifit.h"
      #include "gsl/gsl_multilarge.h"
      #include "gsl/gsl_multifit_nlinear.h"
      #include "gsl/gsl_multilarge_nlinear.h"
    }
  }

  /* Minimization

     This module considers both 1D and multi-D minimization routines
  */
  module Minimize {
    use Common;
    use Array;

    extern {
      #include "gsl/gsl_min.h"
      #include "gsl/gsl_multimin.h"
    }
  }

  /* Monte Carlo Integration

     Includes
       - ``gsl_monte_plain.h``
       - ``gsl_monte_miser.h``
       - ``gsl_monte_vegas.h``
  */
  module MonteCarlo {
    extern {
      #include "gsl/gsl_monte_plain.h"
      #include "gsl/gsl_monte_miser.h"
      #include "gsl/gsl_monte_vegas.h"
    }
  }


  /* Multisets

     Includes ``gsl_multiset.h``
  */
  module Multiset {
    extern {
      #include "gsl/gsl_multiset.h"
    }
  }

  /* N-tuples

     Includes ``gsl_ntuple.h``
  */
  module Ntuple {
    extern {
      #include "gsl/gsl_ntuple.h"
    }
  }


  /* Ordinary Differential Equations

     Include "gsl_odeiv2.h"
  */
  module ODE {
    use Common;

    extern {
      #include "gsl/gsl_odeiv2.h"
    }
  }

  /* Permutations

     Includes ``gsl_permutation.h``
  */
  module Permutations {
    extern {
      #include "gsl/gsl_permutation.h"
    }
  }

  /* Polynomials

     Includes ``gsl_poly.h``
  */
  module Polynomials {
    extern {
      #include "gsl/gsl_poly.h"
    }
  }

  /* Random number generation

     Based on ``gsl_rng.h``
  */
  module RNG {
    extern {
      #include "gsl/gsl_rng.h"
    }

    // Not a complete list as yet??
    enum RNGAlgorithms {
      MT19937,
      MRG
    }

    class Random {
      const rtype : RNGAlgorithms;
      var r : c_ptr(gsl_rng);

      // We make the seed an atomic integer
      // so that, in the future, we can build copy constructors
      // with new types.
      var saveseed : atomic uint;

      proc init(typ : RNGAlgorithms = RNGAlgorithms.MT19937, iseed : uint = 0) {
        rtype = typ;
        this.complete();
        select rtype {
          when RNGAlgorithms.MT19937 do r = gsl_rng_alloc(gsl_rng_mt19937);
          when RNGAlgorithms.MRG do r = gsl_rng_alloc(gsl_rng_mrg);
          otherwise halt("Unknown RNGAlgorithm");
        }
        this.reseed(iseed);
      }

      proc deinit() {
        gsl_rng_free(r);
      }

      proc reseed(iseed : uint) {
        saveseed.write(iseed);
        gsl_rng_set(r, iseed : c_ulong);
      }

      proc seed : uint {
        return saveseed.read();
      }

      proc name : string {
        return gsl_rng_name(r) : string;
      }

      proc getRaw() : c_ulong {
        return gsl_rng_get(r);
      }

      // Return a randomly sampled real number.
      // If pos is true, then the returned value is in (0,1)
      // otherwise it is in [0,1)
      proc get(pos : bool = false) : real {
        if pos {
          return gsl_rng_uniform_pos(r);
        } else {
          return gsl_rng_uniform(r);
        }
      }

      // Return a randomly sampled integer in [0,n)
      proc getInt(n) : uint(64) {
        return gsl_rng_uniform_int(r,n:c_ulong):uint(64);
      }



      // End Class Random
    }



  }

  /* Random number distributions.

     Based on ``gsl_randist.h`` and ``gsl_cdf.h``
  */
  module RanDist {
    use RNG;

    extern {
      #include "gsl/gsl_randist.h"
      #include "gsl/gsl_cdf.h"
    }

    // Wrapper for Gaussian distributions
    record Gaussian {
      var sigma : real;

      proc init(sig=1.0) {
        sigma = sig;
      }

      proc get(engine : Random) : real {
        return gsl_ran_gaussian_ziggurat(engine.r, sigma);
      }

      proc pdf(x) : real {
        return gsl_ran_gaussian_pdf(x : c_double, sigma);
      }

      proc P(x) : real {
        return gsl_cdf_gaussian_P(x : c_double, sigma);
      }

      proc Q(x) : real {
        return gsl_cdf_gaussian_Q(x : c_double, sigma);
      }


      proc Pinv(x) : real {
        return gsl_cdf_gaussian_Pinv(x : c_double, sigma);
      }

      proc Qinv(x) : real {
        return gsl_cdf_gaussian_Qinv(x : c_double, sigma);
      }
    }


    // Wrapper for chi-squared distributions
    record ChiSquared {
      var nu : real;

      proc get(engine : Random) : real {
        return gsl_ran_chisq(engine.r, nu);
      }

      proc pdf(x) : real {
        return gsl_ran_chisq_pdf(x : c_double, nu);
      }

      proc P(x) : real {
        return gsl_cdf_chisq_P(x : c_double, nu);
      }

      proc Q(x) : real {
        return gsl_cdf_chisq_Q(x : c_double, nu);
      }


      proc Pinv(x) : real {
        return gsl_cdf_chisq_Pinv(x : c_double, nu);
      }

      proc Qinv(x) : real {
        return gsl_cdf_chisq_Qinv(x : c_double, nu);
      }
    }

    // End Random module
  }

  /* Root Finding

     Includes ``gsl_roots.h`` and ``gsl_multiroots.h``
  */
  module Roots {
    use Common;
    extern {
      #include "gsl/gsl_roots.h"
      #include "gsl/gsl_multiroots.h"
    }
  }

  /* Quasi-Random number generation

     Based on ``gsl_qrng.h``
  */
  module QRNG {
    extern {
      #include "gsl/gsl_qrng.h"
    }
  }

  /* Series Acceleration

     Include ``gsl_sum.h``
  */
  module SeriesAccel {
    extern {
      #include "gsl/gsl_sum.h"
    }
  }

  /* Simulated Annealing

     Includes ``gsl_siman.h``
  */
  module SimAnneal {
    extern {
      #include "gsl/gsl_siman.h"
    }
  }


  /* Sorting

     Includes ``gsl_heapsort.h``, ``gsl_sort.h`` and
     ``gsl_sort_vector.h``.
   */
  module Sorting {
    use Array;

    extern {
      #include "gsl/gsl_heapsort.h"
      #include "gsl/gsl_sort.h"
      #include "gsl/gsl_sort_vector.h"
    }
  }

  /* Sparse Matrix support

     Includes sparse matrices, BLAS and linear algebra.
  */
  module Sparse {
    use Common;
    use Array;

    extern {
      #include "gsl/gsl_spmatrix.h"
      #include "gsl/gsl_spblas.h"
      #include "gsl/gsl_splinalg.h"
    }
  }


  /* Special functions.

     This includes ``gsl_sf.h`` and provides access
     to all the GSL special functions.
  */
  module SpecFun {
    extern {
      #include "gsl/gsl_sf.h"
    }
  }

  /* Statistics

     This include ``gsl_statistics_double.h``
  */
  module Statistics {
    extern {
      #include "gsl/gsl_statistics_double.h"
      #include "gsl/gsl_rstat.h"
      #include "gsl/gsl_movstat.h"
    }
  }

}