/** \file DefaultVelTraits.hpp
 *  \brief Traits for the default velocity field
 */

#ifndef DEFAULTVELTRAITS_HPP
#define DEFAULTVELTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/DefaultVelFieldTraits.hpp"
#include "PhysicalFields/Domains/FullSphereField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the default velocity field
    */
   class DefaultVelTraits
   {
      public:
         /// Typedef for the variable type
         typedef FullSphereField<DefaultVelFieldTraits> VariableType;
   };

}

#endif // DEFAULTVELTRAITS_HPP
