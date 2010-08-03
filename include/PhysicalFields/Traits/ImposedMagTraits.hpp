/** \file ImposedMagTraits.hpp
 *  \brief Traits for the imposed magnetic field
 */

#ifndef IMPOSEDMAGTRAITS_HPP
#define IMPOSEDMAGTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/ImposedMagFieldTraits.hpp"
#include "PhysicalFields/Domains/FullSphereField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the imposed magnetic field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class ImposedMagTraits
   {
      public:
         /// Typedef for the variable type
         typedef FullSphereField<TSimType, ImposedMagFieldTraits> VariableType;
   };

}

#endif // IMPOSEDMAGTRAITS_HPP
