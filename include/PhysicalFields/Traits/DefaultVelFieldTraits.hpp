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
#include "PhysicalFields/Fields/PhysicalTorPolBase.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolField.hpp"

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
         typedef typename PhysicalTorPolField<TSimType, PhysicalTorPolBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolField<TSimType, PhysicalTorPolBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolField<TSimType, PhysicalTorPolBase> OCFieldType;
   };

}

#endif // DEFAULTVELFIELDTRAITS_HPP
