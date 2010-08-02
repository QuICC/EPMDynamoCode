/** \file DefaultCodTraits.hpp
 *  \brief Traits for the default codensity scalar field
 */

#ifndef DEFAULTCODTRAITS_HPP
#define DEFAULTCODTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/DefaultCodDomainTraits.hpp"
#include "PhysicalFields/Sources/OuterCoreSource.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the default codensity scalar field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class DefaultCodTraits
   {
      public:
         /// Typedef for the domain traits
         typedef DefaultCodDomainTraits<TSimType> DomainTraits;

         /// Typedef for the variable type
         typedef OuterCoreSource<TSimType, DefaultCodDomainTraits> VariableType;
   };

}

#endif // DEFAULTCODTRAITS_HPP
