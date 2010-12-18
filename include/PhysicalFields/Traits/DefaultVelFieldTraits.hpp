/** \file DefaultVelFieldTraits.hpp
 *  \brief Traits for the field of the default velocity field
 */

#ifndef DEFAULTVELFIELDTRAITS_HPP
#define DEFAULTVELFIELDTRAITS_HPP

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
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class DefaultVelFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalTorPolCurl<TSimType, PhysicalTorPolCurlBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolCurl<TSimType, PhysicalTorPolCurlBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolCurl<TSimType, PhysicalTorPolCurlBase> OCFieldType;
   };

}

#endif // DEFAULTVELFIELDTRAITS_HPP
