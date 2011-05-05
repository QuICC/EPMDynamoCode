/** \file FullCoreField.hpp
 *  \brief Implementation for a full core field
 */

#ifndef FULLCOREFIELD_HPP
#define FULLCOREFIELD_HPP

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
    * \brief Implementation for a full core field
    *
    * \tparam TFieldTraits Traits to describe the field
    */
   template <typename TFieldTraits> class FullCoreField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the inner core field type
         typedef typename TFieldTraits::ICFieldType  ICFieldType;

         /// Typedef for the outer core field type
         typedef typename TFieldTraits::OCFieldType  OCFieldType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Transform object
          */
         FullCoreField(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~FullCoreField() {};

         /**
          * @brief Get Physical variable in IC
          */
         const ICFieldType&   ic() const;

         /**
          * @brief Set Physical variable in IC
          */
         ICFieldType&   rIc();

         /**
          * @brief Get Physical variable in OC
          */
         const OCFieldType&   oc() const;

         /**
          * @brief Set Physical variable in OC
          */
         OCFieldType&   rOc();

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
      protected:
         /**
          * @brief Storage for the IC part
          */
         ICFieldType   mICVar;

         /**
          * @brief Storage for the OC part
          */
         OCFieldType   mOCVar;

      private:
   };

   template <typename TFieldTraits> FullCoreField<TFieldTraits>::FullCoreField(SmartTruncation pTrunc, typename FullCoreField<TFieldTraits>::TransformType &transform)
      : mICVar(pTrunc, transform), mOCVar(pTrunc, transform)
   {
   }

   template <typename TFieldTraits> const FullCoreField<TFieldTraits>::ICFieldType& FullCoreField<TFieldTraits>::ic() const
   {
      return this->mICVar;
   }

   template <typename TFieldTraits> FullCoreField<TFieldTraits>::ICFieldType& FullCoreField<TFieldTraits>::rIc()
   {
      return this->mICVar;
   }

   template <typename TFieldTraits> const FullCoreField<TFieldTraits>::OCFieldType& FullCoreField<TFieldTraits>::oc() const
   {
      return this->mOCVar;
   }

   template <typename TFieldTraits> FullCoreField<TFieldTraits>::OCFieldType& FullCoreField<TFieldTraits>::rOc()
   {
      return this->mOCVar;
   }

   template <typename TFieldTraits> void FullCoreField<TFieldTraits>::initialiseZeros()
   {
      this->mICVar.initialiseZeros();

      this->mOCVar.initialiseZeros();
   }

}

#endif // FULLCOREFIELD_HPP
