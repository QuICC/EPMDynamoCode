/** \file GenVelFieldTraits.hpp
 *  \brief Traits for the field of the generators velocity field
 */

#ifndef GENVELFIELDTRAITS_HPP
#define GENVELFIELDTRAITS_HPP

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
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class GenVelFieldTraits
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

#endif // GENVELFIELDTRAITS_HPP
