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
#include "PhysicalFields/Fields/PhysicalTorPolBase.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolField.hpp"
#include "PhysicalFields/Imposed/ImposedTorPol.hpp"

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
         typedef typename PhysicalTorPolField<TSimType, PhysicalTorPolBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolField<TSimType, PhysicalTorPolBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolField<TSimType, ImposedTorPol> OCFieldType;
   };

}

#endif // GENMAGFIELDTRAITS_HPP
