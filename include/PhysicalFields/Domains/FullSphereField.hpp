/** \file FullSphereField.hpp
 *  \brief Implementation for a full sphere field
 */

#ifndef FULLSPHEREFIELD_HPP
#define FULLSPHEREFIELD_HPP

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
    * \brief Implementation for a full sphere field
    *
    * \tparam TFieldTraits Traits to describe the field
    */
   template <typename TFieldTraits> class FullSphereField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the Field type
         typedef typename TFieldTraits::OCFieldType  OCFieldType;

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

   template <typename TFieldTraits> FullSphereField<TFieldTraits>::FullSphereField(SmartTruncation pTrunc, typename FullSphereField<TFieldTraits>::TransformType &transform)
      : mOCVar(pTrunc, transform)
   {
   }

   template <typename TFieldTraits> const typename FullSphereField<TFieldTraits>::OCFieldType& FullSphereField<TFieldTraits>::oc() const
   {
      return this->mOCVar;
   }

   template <typename TFieldTraits> typename FullSphereField<TFieldTraits>::OCFieldType& FullSphereField<TFieldTraits>::rOc()
   {
      return this->mOCVar;
   }

   template <typename TFieldTraits> void  FullSphereField<TFieldTraits>::initialiseZeros()
   {
      this->mOCVar.initialiseZeros();
   }

}

#endif // FULLSPHEREFIELD_HPP
