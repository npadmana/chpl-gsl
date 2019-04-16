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