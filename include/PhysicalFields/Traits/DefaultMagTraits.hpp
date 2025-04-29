/** \file DefaultMagTraits.hpp
 *  \brief Traits for the default magnetic field
 */

#ifndef DEFAULTMAGTRAITS_HPP
#define DEFAULTMAGTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/DefaultMagFieldTraits.hpp"
#include "PhysicalFields/Domains/FullSphereField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the default magnetic field
    */
   class DefaultMagTraits
   {
      public:
         /// Typedef for the variable type
         typedef FullSphereField<DefaultMagFieldTraits> VariableType;
   };

}

#endif // DEFAULTMAGTRAITS_HPP
