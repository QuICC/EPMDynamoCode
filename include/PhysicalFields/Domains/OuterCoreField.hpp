/** \file OuterCoreField.hpp
 *  \brief Implementation for an outer core field
 */

#ifndef OUTERCOREFIELD_HPP
#define OUTERCOREFIELD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation for an outer core field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldTraits Type of the field
    */
   template <typename TSimType, template <typename> class TFieldTraits> class OuterCoreField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the outer core field type
         typedef typename TFieldTraits<TSimType>::OCFieldType  OCFieldType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Transform object
          */
         OuterCoreField(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~OuterCoreField() {};

         /**
          * @brief Get Physical variable in OC (or full sphere)
          */
         const OCFieldType&   oc() const;

         /**
          * @brief Set Physical variable in OC (or full sphere)
          */
         OCFieldType&   rOc();
         
      protected:
         /**
          * @brief Storage for the OC (or full sphere) part
          */
         OCFieldType   mOCVar;

      private:
   };

   template <typename TSimType, template <typename>  class TFieldTraits> OuterCoreField<TSimType, TFieldTraits>::OuterCoreField(SmartTruncation pTrunc, typename OuterCoreField<TSimType,TFieldTraits>::TransformType &transform)
      : mOCVar(pTrunc, transform)
   {
   }

   template <typename TSimType, template <typename>  class TFieldTraits> const OuterCoreField<TSimType, TFieldTraits>::OCFieldType& OuterCoreField<TSimType, TFieldTraits>::oc() const
   {
      return this->mOCVar;
   }

   template <typename TSimType, template <typename>  class TFieldTraits> OuterCoreField<TSimType, TFieldTraits>::OCFieldType& OuterCoreField<TSimType, TFieldTraits>::rOc()
   {
      return this->mOCVar;
   }

}

#endif // OUTERCOREFIELD_HPP
