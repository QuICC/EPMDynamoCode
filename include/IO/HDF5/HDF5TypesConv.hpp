/** \file HDF5TypesConv.hpp
 *  \brief Implementation of a simple HDF5 datatype converter
 */

#ifndef HDF5TYPESCONV_HPP
#define HDF5TYPESCONV_HPP

// System includes
//
#include <complex>
#include <hdf5.h>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief This class implements a few type conversion routines for the HDF5 files
    */
   class HDF5TypesConv
   {
      public:

         /**
          * @brief Simple empty constructor
          */
         HDF5TypesConv();

         /**
          * @brief Simple empty destructor
          */
         virtual ~HDF5TypesConv() {};

         /**
          * @brief Get HDF5 predefined type for a scalar type
          *
          * \tparam T type to convert to HDF5 type
          */
         template <typename T> hid_t ptype();

         /**
          * @brief Get HDF5 predefined type for an array type
          *
          * \tparam T type to convert to HDF5 type
          */
         template <typename T> hid_t atype();
   };

   template <> inline hid_t HDF5TypesConv::ptype<int>()
   {
      return H5T_NATIVE_INT;
   }

   template <> inline hid_t HDF5TypesConv::ptype<float>()
   {
      return H5T_NATIVE_FLOAT;
   }

   template <> inline hid_t HDF5TypesConv::ptype<double>()
   {
      return H5T_NATIVE_DOUBLE;
   }

   /// Specialised method for a complex<float>
   template <> inline hid_t HDF5TypesConv::atype<std::complex<float> >()
   {
      hsize_t dims = 2;
      return H5Tarray_create(H5T_NATIVE_FLOAT, 1, &dims);
   }

   /// Specialised method for a complex<double>
   template <> inline hid_t HDF5TypesConv::atype<std::complex<double> >()
   {
      hsize_t dims = 2;
      return H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &dims);
   }

}

#endif // HDF5TYPESCONV_HPP
