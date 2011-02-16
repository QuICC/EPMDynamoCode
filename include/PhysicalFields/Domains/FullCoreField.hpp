/** \file FullCoreField.hpp
 *  \brief Implementation for a full core field
 */

#ifndef FULLCOREFIELD_HPP
#define FULLCOREFIELD_HPP

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
    * \brief Implementation for a full core field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldTraits Traits to describe the field
    */
   template <typename TSimType, template <typename> class TFieldTraits> class FullCoreField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the inner core field type
         typedef typename TFieldTraits<TSimType>::ICFieldType  ICFieldType;

         /// Typedef for the outer core field type
         typedef typename TFieldTraits<TSimType>::OCFieldType  OCFieldType;

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

   template <typename TSimType, template <typename>  class TFieldTraits> FullCoreField<TSimType, TFieldTraits>::FullCoreField(SmartTruncation pTrunc, typename FullCoreField<TSimType,TFieldTraits>::TransformType &transform)
      : mICVar(pTrunc, transform), mOCVar(pTrunc, transform)
   {
   }

   template <typename TSimType, template <typename>  class TFieldTraits> const FullCoreField<TSimType, TFieldTraits>::ICFieldType& FullCoreField<TSimType, TFieldTraits>::ic() const
   {
      return this->mICVar;
   }

   template <typename TSimType, template <typename>  class TFieldTraits> FullCoreField<TSimType, TFieldTraits>::ICFieldType& FullCoreField<TSimType, TFieldTraits>::rIc()
   {
      return this->mICVar;
   }

   template <typename TSimType, template <typename>  class TFieldTraits> const FullCoreField<TSimType, TFieldTraits>::OCFieldType& FullCoreField<TSimType, TFieldTraits>::oc() const
   {
      return this->mOCVar;
   }

   template <typename TSimType, template <typename>  class TFieldTraits> FullCoreField<TSimType, TFieldTraits>::OCFieldType& FullCoreField<TSimType, TFieldTraits>::rOc()
   {
      return this->mOCVar;
   }

   template <typename TSimType, template <typename>  class TFieldTraits> void FullCoreField<TSimType, TFieldTraits>::initialiseZeros()
   {
      this->mICVar.initialiseZeros();

      this->mOCVar.initialiseZeros();
   }

}

#endif // FULLCOREFIELD_HPP
