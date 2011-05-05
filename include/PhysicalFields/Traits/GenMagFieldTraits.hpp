/** \file GenMagFieldTraits.hpp
 *  \brief Traits for the field of the generators magnetic field
 */

#ifndef GENMAGFIELDTRAITS_HPP
#define GENMAGFIELDTRAITS_HPP

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
    * \brief Traits for the field of the generators magnetic field
    */
   class GenMagFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalTorPolCurlField<PhysicalTorPolCurlBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolCurl<PhysicalTorPolCurlBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolCurl<ImposedTorPol> OCFieldType;
   };

}

#endif // GENMAGFIELDTRAITS_HPP
