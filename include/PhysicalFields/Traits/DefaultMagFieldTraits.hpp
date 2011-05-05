/** \file DefaultMagFieldTraits.hpp
 *  \brief Traits for the field of the default magnetic field
 */

#ifndef DEFAULTMAGFIELDTRAITS_HPP
#define DEFAULTMAGFIELDTRAITS_HPP

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
    * \brief Traits for the field of the default magnetic field
    */
   class DefaultMagFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalTorPolCurl<PhysicalTorPolCurlBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolCurl<PhysicalTorPolCurlBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolCurl<PhysicalTorPolCurlBase> OCFieldType;
   };

}

#endif // DEFAULTMAGFIELDTRAITS_HPP
