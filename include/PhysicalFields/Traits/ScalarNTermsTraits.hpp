/** \file ScalarNTermsTraits.hpp
 *  \brief Traits for the default scalar non linear terms
 */

#ifndef SCALARNTERMSTRAITS_HPP
#define SCALARNTERMSTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/ScalarNTermsFieldTraits.hpp"
#include "PhysicalFields/Domains/FullSphereField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the default scalar non linear terms
    */
   class ScalarNTermsTraits
   {
      public:
         /// Typedef for the variable type
         typedef FullSphereField<ScalarNTermsFieldTraits> VariableType;
   };

}

#endif // SCALARNTERMSTRAITS_HPP
