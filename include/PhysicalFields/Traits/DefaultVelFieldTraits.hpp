/** \file DefaultVelFieldTraits.hpp
 *  \brief Traits for the field of the default velocity field
 */

#ifndef DEFAULTVELFIELDTRAITS_HPP
#define DEFAULTVELFIELDTRAITS_HPP

// Configuraton includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalTorPolCurlBase.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolCurl.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the field of the default velocity field
    */
   class DefaultVelFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef PhysicalTorPolCurl<PhysicalTorPolCurlBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolCurl<PhysicalTorPolCurlBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolCurl<PhysicalTorPolCurlBase> OCFieldType;
   };

}

#endif // DEFAULTVELFIELDTRAITS_HPP
