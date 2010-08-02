/** \file DefaultVelTraits.hpp
 *  \brief Traits for the default velocity field
 */

#ifndef DEFAULTVELTRAITS_HPP
#define DEFAULTVELTRAITS_HPP

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
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class DefaultVelTraits
   {
      public:
         /// Typedef for the magnetic fields
         typedef  DefaultVelFieldTraits<TSimType>   FieldTraits;

         /// Typedef for the variable type
         typedef FullSphereField<TSimType, DefaultVelFieldTraits> VariableType;
   };

}

#endif // DEFAULTVELTRAITS_HPP
