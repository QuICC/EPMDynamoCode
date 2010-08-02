/** \file DefaultCodFieldTraits.hpp
 *  \brief Traits for the fields of the default codensitfy scalar field
 */

#ifndef DEFAULTCODFIELDTRAITS_HPP
#define DEFAULTCODFIELDTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalScalarBase.hpp"
#include "PhysicalFields/Fields/PhysicalScalar.hpp"

namespace EPMDynamo {

   /**
    *  \brief Traits for the fields of the default codensitfy scalar field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class DefaultCodFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef typename PhysicalScalar<TSimType, PhysicalScalarBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the inner core field
         typedef PhysicalScalar<TSimType, PhysicalScalarBase> ICFieldType;

         /// Typedef for the outer core field
         typedef PhysicalScalar<TSimType, PhysicalScalarBase> OCFieldType;
   };

}

#endif // DEFAULTCODFIELDTRAITS_HPP
