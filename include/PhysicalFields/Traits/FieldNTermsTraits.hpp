/** \file FieldNTermsTraits.hpp
 *  \brief Traits for the vector field non linear terms
 */

#ifndef FIELDNTERMSTRAITS_HPP
#define FIELDNTERMSTRAITS_HPP

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
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType>  class FieldNTermsTraits
   {
      public:
         /// Typedef for the variable type
         typedef FullSphereField<TSimType, FieldNTermsFieldTraits> VariableType;
   };

}

#endif // FIELDNTERMSTRAITS_HPP
