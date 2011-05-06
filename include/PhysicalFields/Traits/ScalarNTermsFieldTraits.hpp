/** \file ScalarNTermsFieldTraits.hpp
 *  \brief Traits for the fields of the scalar non linear terms
 */

#ifndef SCALARNTERMSFIELDTRAITS_HPP
#define SCALARNTERMSFIELDTRAITS_HPP

// Configuration includes
//

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
    */
   class ScalarNTermsFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef PhysicalScalar<PhysicalScalarBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the inner core field
         typedef PhysicalScalar<PhysicalScalarBase> ICFieldType;

         /// Typedef for the outer core field
         typedef PhysicalScalar<PhysicalScalarBase> OCFieldType;
   };

}

#endif // SCALARNTERMSFIELDTRAITS_HPP
