/** \file OuterCoreField.hpp
 *  \brief Implementation for an outer core field
 */

#ifndef OUTERCOREFIELD_HPP
#define OUTERCOREFIELD_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation for an outer core field
    *
    * \tparam TFieldTraits Type of the field
    */
   template <typename TFieldTraits> class OuterCoreField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the outer core field type
         typedef typename TFieldTraits::OCFieldType  OCFieldType;

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

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
      protected:
         /**
          * @brief Storage for the OC (or full sphere) part
          */
         OCFieldType   mOCVar;

      private:
   };

   template <typename TFieldTraits> OuterCoreField<TFieldTraits>::OuterCoreField(SmartTruncation pTrunc, typename OuterCoreField<TFieldTraits>::TransformType &transform)
      : mOCVar(pTrunc, transform)
   {
   }

   template <typename TFieldTraits> const OuterCoreField<TFieldTraits>::OCFieldType& OuterCoreField<TFieldTraits>::oc() const
   {
      return this->mOCVar;
   }

   template <typename TFieldTraits> OuterCoreField<TFieldTraits>::OCFieldType& OuterCoreField<TFieldTraits>::rOc()
   {
      return this->mOCVar;
   }

   template <typename TFieldTraits>  OuterCoreField<TFieldTraits>::initialiseZeros()
   {
      this->mOCVar.initialiseZeros();
   }

}

#endif // OUTERCOREFIELD_HPP
