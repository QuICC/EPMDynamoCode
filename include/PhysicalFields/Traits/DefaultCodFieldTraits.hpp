/** \file DefaultCodFieldTraits.hpp
 *  \brief Traits for the fields of the default codensitfy scalar field
 */

#ifndef DEFAULTCODFIELDTRAITS_HPP
#define DEFAULTCODFIELDTRAITS_HPP

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

namespace EPMDynamo {

   /**
    *  \brief Traits for the fields of the default codensitfy scalar field
    */
   class DefaultCodFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef PhysicalScalarGradient<PhysicalScalarGradientBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the inner core field
         typedef PhysicalScalarGradient<PhysicalScalarGradientBase> ICFieldType;

         /// Typedef for the outer core field
         typedef PhysicalScalarGradient<PhysicalScalarGradientBase> OCFieldType;
   };

}

#endif // DEFAULTCODFIELDTRAITS_HPP
