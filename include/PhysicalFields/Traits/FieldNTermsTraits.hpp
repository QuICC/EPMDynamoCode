/** \file FieldNTermsTraits.hpp
 *  \brief Traits for the vector field non linear terms
 */

#ifndef FIELDNTERMSTRAITS_HPP
#define FIELDNTERMSTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Traits/FieldNTermsFieldTraits.hpp"
#include "PhysicalFields/Domains/FullSphereField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the vector field non linear terms
    */
   class FieldNTermsTraits
   {
      public:
         /// Typedef for the variable type
         typedef FullSphereField<FieldNTermsFieldTraits> VariableType;
   };

}

#endif // FIELDNTERMSTRAITS_HPP
