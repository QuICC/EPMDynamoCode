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
#include "PhysicalFields/Fields/PhysicalTorPolBase.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolField.hpp"
#include "PhysicalFields/Imposed/ImposedTorPol.hpp"

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
         typedef typename PhysicalTorPolField<TSimType, PhysicalTorPolBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolField<TSimType, PhysicalTorPolBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolField<TSimType, ImposedTorPol> OCFieldType;
   };

}

#endif // GENVELFIELDTRAITS_HPP
