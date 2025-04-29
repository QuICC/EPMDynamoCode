/** \file InnerCoreField.hpp
 *  \brief Implementation for an inner core field
 */

#ifndef INNERCOREFIELD_HPP
#define INNERCOREFIELD_HPP

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
    * \brief Implementation for an inner core field
    *
    * \tparam TFieldTraits Type of the field
    */
   template <typename TFieldTraits> class InnerCoreField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the inner core field type
         typedef typename TFieldTraits::ICFieldType  ICFieldType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Transform object
          */
         InnerCoreField(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~InnerCoreField() {};

         /**
          * @brief Get Physical variable in IC
          */
         const ICFieldType&   ic() const;

         /**
          * @brief Set Physical variable in IC
          */
         ICFieldType&   rIc();

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
      protected:
         /**
          * @brief Storage for the IC part
          */
         ICFieldType   mICVar;

      private:
   };

   template <typename TFieldTraits> InnerCoreField<TFieldTraits>::InnerCoreField(SmartTruncation pTrunc, typename InnerCoreField<TFieldTraits>::TransformType &transform)
      : mICVar(pTrunc, transform)
   {
   }

   template <typename TFieldTraits> const InnerCoreField<TFieldTraits>::ICFieldType& InnerCoreField<TFieldTraits>::ic() const
   {
      return this->mICVar;
   }

   template <typename TFieldTraits> InnerCoreField<TFieldTraits>::ICFieldType& InnerCoreField<TFieldTraits>::rIc()
   {
      return this->mICVar;
   }

   template <typename TFieldTraits> void InnerCoreField<TFieldTraits>::initialiseZeros()
   {
      this->mICVar.initialiseZeros();
   }

}

#endif // INNERCOREFIELD_HPP
