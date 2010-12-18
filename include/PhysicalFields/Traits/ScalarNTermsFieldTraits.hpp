/** \file ScalarNTermsFieldTraits.hpp
 *  \brief Traits for the fields of the scalar non linear terms
 */

#ifndef SCALARNTERMSFIELDTRAITS_HPP
#define SCALARNTERMSFIELDTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalScalarBase.hpp"
#include "PhysicalFields/Fields/PhysicalScalar.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the fields of the scalar non linear terms
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class ScalarNTermsFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalScalar<TSimType, PhysicalScalarBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the inner core field
         typedef PhysicalScalar<TSimType, PhysicalScalarBase> ICFieldType;

         /// Typedef for the outer core field
         typedef PhysicalScalar<TSimType, PhysicalScalarBase> OCFieldType;
   };

}

#endif // SCALARNTERMSFIELDTRAITS_HPP
