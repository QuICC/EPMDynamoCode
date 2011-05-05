/** \file ImposedMagTraits.hpp
 *  \brief Traits for the imposed magnetic field
 */

#ifndef IMPOSEDMAGTRAITS_HPP
#define IMPOSEDMAGTRAITS_HPP

// Configuration includes
//

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
    */
   class ImposedMagTraits
   {
      public:
         /// Typedef for the variable type
         typedef FullSphereField<ImposedMagFieldTraits> VariableType;
   };

}

#endif // IMPOSEDMAGTRAITS_HPP
