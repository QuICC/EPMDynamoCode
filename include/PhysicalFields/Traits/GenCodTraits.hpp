/** \file GenCodTraits.hpp
 *  \brief Traits for the generators codensity scalar field
 */

#ifndef GENCODTRAITS_HPP
#define GENCODTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/GenCodDomainTraits.hpp"
#include "PhysicalFields/Sources/OuterCoreSource.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the generators codensity scalar field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class GenCodTraits
   {
      public:
         /// Typedef for the variable type
         typedef OuterCoreSource<TSimType, GenCodDomainTraits> VariableType;
   };

}

#endif // GENCODTRAITS_HPP
