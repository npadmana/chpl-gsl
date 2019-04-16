/*
  This is the top-level GSL module. Most of the GSL functionality
  is included in a series of submodules, which are not automatically
  included.

  Instead of explicit declarations, this makes heavy use of extern blocks. The
  Chapel code here is restricted to convenience functions and workarounds.

  The mapping from sub-modules to GSL chapters is :

  * :mod:`Complex` : Complex numbers
  * :mod:`RNG` : Random number generation
  * :mod:`SpecFun` : Special Functions
  
*/
module GSL {

  use Common;

  module Common {

    extern record gsl_function {
      var function : c_fn_ptr;
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
      var p : c_ptr(gsl_matrix); // Pointer to GSL Matrix
      var pdata : c_ptr(c_double); // Pointer to the actual data
      var D : domain(2);

      proc init(m : size_t, n : size_t) {
        p = gsl_matrix_alloc(m, n);
        pdata = (p.deref()).data;
        D = {0.. #m, 0.. #n};
      }

      proc deinit() {
        gsl_matrix_free(p);
      }

      proc access(i : size_t, j : size_t) ref {
        return (gsl_matrix_ptr(p, i, j)).deref();
      }

      proc this(i : int, j : int) ref {
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
      type gslVecType;
      var p : gslVecType; // Pointer to GSL Vector
      var pdata : c_ptr(eltType); // Pointer to the actual data
      var D : domain(1);

      proc init(n) where isIntegral(n) {
        eltType = c_double;
        gslVecType = c_ptr(gsl_vector);
        this.complete();
        p = gsl_vector_alloc(n:size_t);
        pdata = (p.deref()).data;
        D = {0.. #n};
      }

      proc init(type t, n) where isIntegral(n) {
        eltType = t;
        select t {
            when real(64) do gslVecType = c_ptr(gsl_vector);
            when c_double do gslVecType = c_ptr(gsl_vector);
            when c_int do gslVecType = c_ptr(gsl_vector_int);
            otherwise {
              gslVecType = bool;
              compilerError("Unimplemented type for GSLVector : ",t:string);
            }
          }
        this.complete();
        select t {
            when real(64) do p = gsl_vector_alloc(n:size_t);
            when c_double do p = gsl_vector_alloc(n:size_t);
            when c_int do p = gsl_vector_int_alloc(n:size_t);
          }
        pdata = (p.deref()).data;
        D = {0.. #n};
      }

      proc deinit() {
        select eltType {
            when real(64) do gsl_vector_free(p);
            when c_double do gsl_vector_free(p);
            when c_int do gsl_vector_int_free(p);
          }
      }

      proc access(i : size_t) ref {
        select eltType {
            when real(64) do return (gsl_vector_ptr(p,i)).deref();
            when c_double do return (gsl_vector_ptr(p,i)).deref();
            when c_int do return (gsl_vector_int_ptr(p,i)).deref();
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


  /* Numerical Differentiation.

     Based on ``gsl_deriv.h``.
  */
  module Deriv {
    use Common;

    extern {
      #include "gsl/gsl_deriv.h"
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


  /* Random number generation

     Based on ``gsl_rng.h``
  */
  module RNG {
    extern {
      #include "gsl/gsl_rng.h"
    }
  }

  /* Random number distributions.

     Based on ``gsl_randist.h`` and ``gsl_cdf.h``
  */
  module RanDist {
    extern {
      #include "gsl/gsl_randist.h"
      #include "gsl/gsl_cdf.h"
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

  /* Special functions.

     This includes ``gsl_sf.h`` and provides access
     to all the GSL special functions.
  */
  module SpecFun {
    extern {
      #include "gsl/gsl_sf.h"
    }
  }

}