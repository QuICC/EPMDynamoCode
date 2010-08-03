/** \file GenCodFieldTraits.hpp
 *  \brief Traits for the fields of the generators codensitfy scalar field
 */

#ifndef GENCODFIELDTRAITS_HPP
#define GENCODFIELDTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalScalarBase.hpp"
#include "PhysicalFields/Fields/PhysicalScalar.hpp"
#include "PhysicalFields/Imposed/ImposedScalar.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the fields of the generators codensitfy scalar field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class GenCodFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalScalar<TSimType, PhysicalScalarBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the inner core field
         typedef PhysicalScalar<TSimType, PhysicalScalarBase> ICFieldType;

         /// Typedef for the outer core field
         typedef PhysicalScalar<TSimType, ImposedScalar> OCFieldType;
   };

}

#endif // GENCODFIELDTRAITS_HPP
