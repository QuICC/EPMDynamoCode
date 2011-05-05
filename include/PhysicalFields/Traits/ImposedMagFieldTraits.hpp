/** \file ImposedMagFieldTraits.hpp
 *  \brief Traits for the field of the imposed magnetic field
 */

#ifndef IMPOSEDMAGFIELDTRAITS_HPP
#define IMPOSEDMAGFIELDTRAITS_HPP

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
    * \brief Traits for the field of the imposed magnetic field
    */
   class ImposedMagFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalTorPolCurl<PhysicalTorPolCurlBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolCurl<PhysicalTorPolCurlBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolCurl<ImposedTorPolCurl> OCFieldType;
   };

}

#endif // IMPOSEDMAGFIELDTRAITS_HPP
