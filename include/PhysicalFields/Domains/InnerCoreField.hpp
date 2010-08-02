/** \file InnerCoreField.hpp
 *  \brief Implementation for an inner core field
 */

#ifndef INNERCOREFIELD_HPP
#define INNERCOREFIELD_HPP

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
    * \brief Implementation for an inner core field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldTraits Type of the field
    */
   template <typename TSimType, template <typename> class TFieldTraits> class InnerCoreField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the inner core field type
         typedef typename TFieldTraits<TSimType>::ICFieldType  ICFieldType;

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
         
      protected:
         /**
          * @brief Storage for the IC part
          */
         ICFieldType   mICVar;

      private:
   };

   template <typename TSimType, template <typename>  class TFieldTraits> InnerCoreField<TSimType, TFieldTraits>::InnerCoreField(SmartTruncation pTrunc, typename InnerCoreField<TSimType,TFieldTraits>::TransformType &transform)
      : mICVar(pTrunc, transform)
   {
   }

   template <typename TSimType, template <typename>  class TFieldTraits> const InnerCoreField<TSimType, TFieldTraits>::ICFieldType<TSimType>& InnerCoreField<TSimType, TFieldTraits>::ic() const
   {
      return this->mICVar;
   }

   template <typename TSimType, template <typename>  class TFieldTraits> InnerCoreField<TSimType, TFieldTraits>::ICFieldType<TSimType>& InnerCoreField<TSimType, TFieldTraits>::rIc()
   {
      return this->mICVar;
   }

}

#endif // INNERCOREFIELD_HPP
