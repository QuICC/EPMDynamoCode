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
#include "PhysicalFields/Fields/PhysicalTorPolBase.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolField.hpp"
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
         typedef typename PhysicalTorPolField<TSimType, PhysicalTorPolBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolField<TSimType, PhysicalTorPolBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolField<TSimType, ImposedTorPol> OCFieldType;
   };

}

#endif // IMPOSEDMAGFIELDTRAITS_HPP
