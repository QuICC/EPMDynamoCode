/** \file GenCodFieldTraits.hpp
 *  \brief Traits for the fields of the generators codensitfy scalar field
 */

#ifndef GENCODFIELDTRAITS_HPP
#define GENCODFIELDTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalScalarGradientBase.hpp"
#include "PhysicalFields/Fields/PhysicalScalarGradient.hpp"
#include "PhysicalFields/Imposed/ImposedScalar.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the fields of the generators codensitfy scalar field
    */
   class GenCodFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalScalarGradient<PhysicalScalarGradientBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the inner core field
         typedef PhysicalScalarGradient<PhysicalScalarGradientBase> ICFieldType;

         /// Typedef for the outer core field
         typedef PhysicalScalarGradient<ImposedScalar> OCFieldType;
   };

}

#endif // GENCODFIELDTRAITS_HPP
