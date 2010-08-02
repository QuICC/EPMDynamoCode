/** \file FullSphereField.hpp
 *  \brief Implementation for a full sphere field
 */

#ifndef FULLSPHEREFIELD_HPP
#define FULLSPHEREFIELD_HPP

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
    * \brief Implementation for a full sphere field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldTraits Traits to describe the field
    */
   template <typename TSimType, template <typename> class TFieldTraits> class FullSphereField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the Field type
         typedef typename TFieldTraits<TSimType>::OCFieldType  OCFieldType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Transform object
          */
         FullSphereField(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~FullSphereField() {};

         /**
          * @brief Get Physical variable in full sphere
          */
         const OCFieldType&   oc() const;

         /**
          * @brief Set Physical variable in full sphere
          */
         OCFieldType&   rOc();
         
      protected:
         /**
          * @brief Storage for the OC (or full sphere) part
          */
         OCFieldType   mOCVar;

      private:
   };

   template <typename TSimType, template <typename> class TFieldTraits> FullSphereField<TSimType, TFieldTraits>::FullSphereField(SmartTruncation pTrunc, typename FullSphereField<TSimType, TFieldTraits>::TransformType &transform)
      : mOCVar(pTrunc, transform)
   {
   }

   template <typename TSimType, template <typename> class TFieldTraits> const typename FullSphereField<TSimType, TFieldTraits>::OCFieldType& FullSphereField<TSimType, TFieldTraits>::oc() const
   {
      return this->mOCVar;
   }

   template <typename TSimType, template <typename> class TFieldTraits> typename FullSphereField<TSimType, TFieldTraits>::OCFieldType& FullSphereField<TSimType, TFieldTraits>::rOc()
   {
      return this->mOCVar;
   }

}

#endif // FULLSPHEREFIELD_HPP
