/** \file DefaultMagTraits.hpp
 *  \brief Traits for the default magnetic field
 */

#ifndef DEFAULTMAGTRAITS_HPP
#define DEFAULTMAGTRAITS_HPP

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
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class DefaultMagTraits
   {
      public:
         /// Typedef for the variable type
         typedef FullSphereField<TSimType, DefaultMagFieldTraits> VariableType;
   };

}

#endif // DEFAULTMAGTRAITS_HPP
