/** \file GenMagFieldTraits.hpp
 *  \brief Traits for the field of the generators magnetic field
 */

#ifndef GENMAGFIELDTRAITS_HPP
#define GENMAGFIELDTRAITS_HPP

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
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class GenMagFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalTorPolCurlField<TSimType, PhysicalTorPolCurlBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolCurl<TSimType, PhysicalTorPolCurlBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolCurl<TSimType, ImposedTorPol> OCFieldType;
   };

}

#endif // GENMAGFIELDTRAITS_HPP
