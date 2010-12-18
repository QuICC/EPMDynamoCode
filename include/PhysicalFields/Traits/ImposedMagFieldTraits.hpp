/** \file ImposedMagFieldTraits.hpp
 *  \brief Traits for the field of the imposed magnetic field
 */

#ifndef IMPOSEDMAGFIELDTRAITS_HPP
#define IMPOSEDMAGFIELDTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalTorPolCurlBase.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolCurl.hpp"
#include "PhysicalFields/Imposed/ImposedTorPol.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the field of the imposed magnetic field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class ImposedMagFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalTorPolCurl<TSimType, PhysicalTorPolCurlBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolCurl<TSimType, PhysicalTorPolCurlBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolCurl<TSimType, ImposedTorPol> OCFieldType;
   };

}

#endif // IMPOSEDMAGFIELDTRAITS_HPP
