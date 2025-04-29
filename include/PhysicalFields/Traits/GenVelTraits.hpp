/** \file GenVelTraits.hpp
 *  \brief Traits for the generators velocity field
 */

#ifndef GENVELTRAITS_HPP
#define GENVELTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/GenVelDomainTraits.hpp"
#include "PhysicalFields/Sources/OuterCoreSource.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the generators velocity field
    */
   class GenVelTraits
   {
      public:
         /// Typedef for the variable type
         typedef OuterCoreSource<GenVelDomainTraits> VariableType;
   };

}

#endif // GENVELTRAITS_HPP
