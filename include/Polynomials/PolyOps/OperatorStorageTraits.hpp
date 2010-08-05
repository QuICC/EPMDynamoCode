/** \file OperatorStorageTraits.hpp
 */

#ifndef OPERATORSTORAGETRAITS_HPP
#define OPERATORSTORAGETRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Polynomials/PolyOps/StoredOperator.hpp"
#include "Polynomials/PolyOps/ComputedOperator.hpp"
#include "Polynomials/PolyOps/POperator.hpp"
#include "Polynomials/PolyOps/CPOperator.hpp"

namespace EPMDynamo {

   /**
    * @brief These traits defines the storage type of the operators
    */
   template <typename TExp> class OperatorStorageTraits
   {
      public:
         #ifdef EPMDYNAMO_BIGMEM
            /// Typedef for the used storage type to be stored operators
            typedef StoredOperator<TExp> StorageType;
         #else
            /// Typedef for the used storage type to be computed operators
            typedef ComputedOperator<TExp> StorageType;
         #endif // EPMDYNAMO_BIGMEM
   };

   /**
    * @brief These traits defines the storage type of the operators (Specialized for POperator)
    */
   template <> class OperatorStorageTraits<POperator>
   {
      public:
         /// Typedef for the used storage type
         typedef ComputedOperator<POperator> StorageType;
   };

   /**
    * @brief These traits defines the storage type of the operators (Specialized for DPOperator)
    */
   template <> class OperatorStorageTraits<CPOperator<EPMFloat> >
   {
      public:
         /// Typedef for the used storage type
         typedef ComputedOperator<CPOperator<EPMFloat> > StorageType;
   };

   /**
    * @brief These traits defines the storage type of the operators (Specialized for ZPOperator)
    */
   template <> class OperatorStorageTraits<CPOperator<EPMComplex> >
   {
      public:
         /// Typedef for the used storage type
         typedef ComputedOperator<CPOperator<EPMComplex> > StorageType;
   };

}

#endif // OPERATORSTORAGETRAITS_HPP
