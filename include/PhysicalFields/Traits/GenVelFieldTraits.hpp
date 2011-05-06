/** \file GenVelFieldTraits.hpp
 *  \brief Traits for the field of the generators velocity field
 */

#ifndef GENVELFIELDTRAITS_HPP
#define GENVELFIELDTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalTorPolCurlBase.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolCurl.hpp"
#include "PhysicalFields/Imposed/ImposedTorPolCurl.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the field of the generators velocity field
    */
   class GenVelFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef PhysicalTorPolCurl<PhysicalTorPolCurlBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolCurl<PhysicalTorPolCurlBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolCurl<ImposedTorPol> OCFieldType;
   };

}

#endif // GENVELFIELDTRAITS_HPP
