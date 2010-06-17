/** \file SimpleField.hpp
 *  \brief Implementation for a simple field
 */

#ifndef SIMPLEFIELD_HPP
#define SIMPLEFIELD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation for a simple field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    */
   template <typename TSimType, template <typename> class TFieldType> class SimpleField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Transform object
          */
         SimpleField(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~SimpleField() {};

         /**
          * @brief Get Physical variable in OC (or full sphere)
          */
         const TFieldType<TSimType>&   oc() const;

         /**
          * @brief Set Physical variable in OC (or full sphere)
          */
         TFieldType<TSimType>&   rOc();
         
      protected:
         /**
          * @brief Storage for the OC (or full sphere) part
          */
         TFieldType<TSimType>   mOCVar;

      private:
   };

   template <typename TSimType, template <typename>  class TFieldType> SimpleField<TSimType, TFieldType>::SimpleField(SmartTruncation pTrunc, typename SimpleField<TSimType,TFieldType>::TransformType &transform)
      : mOCVar(pTrunc, transform)
   {
   }

   template <typename TSimType, template <typename>  class TFieldType> const TFieldType<TSimType>& SimpleField<TSimType, TFieldType>::oc() const
   {
      return this->mOCVar;
   }

   template <typename TSimType, template <typename>  class TFieldType> TFieldType<TSimType>& SimpleField<TSimType, TFieldType>::rOc()
   {
      return this->mOCVar;
   }

}

#endif // SIMPLEFIELD_HPP
