/** \file DefaultCodTraits.hpp
 *  \brief Traits for the default codensity scalar field
 */

#ifndef DEFAULTCODTRAITS_HPP
#define DEFAULTCODTRAITS_HPP

// Configuration includes
//

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
    */
   class DefaultCodTraits
   {
      public:
         /// Typedef for the variable type
         typedef OuterCoreSource<DefaultCodDomainTraits> VariableType;
   };

}

#endif // DEFAULTCODTRAITS_HPP
