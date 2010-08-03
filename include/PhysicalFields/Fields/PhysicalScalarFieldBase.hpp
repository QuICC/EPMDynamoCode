/** \file PhysicalScalarFieldBase.hpp
 *  \brief Base of the implementation of a physical scalar
 */

#ifndef PHYSICALSCALARFIELDBASE_HPP
#define PHYSICALSCALARFIELDBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalBase.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/RTPScalar.hpp"
#include "GeneralFields/RTPField.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a physical scalar
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalScalarFieldBase : public PhysicalBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /**
         * @brief Constructs the basic information required for a physical scalar
         *
         * @param pTrunc Truncation information
         * @param transform Transform object
         */
         PhysicalScalarFieldBase(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalScalarFieldBase() {};

         /**
          * @brief Get the RTP scalar values
          */
         const RTPScalar&   rtp() const;

         /**
          * @brief Set the RTP scalar values
          */
         RTPScalar&   rRTP();

         /**
          * @brief Get the RTP gradient field values
          */
         const RTPField&   grad() const;

         /**
          * @brief Compute RTP values of the gradient of the scalar
          *
          * @param step Current step in multistep transform
          */
         virtual void gradTransform(const int step) = 0;
         
      protected:

         /**
          * @brief Set the RTP gradient field values
          */
         RTPField&   rGrad();

      private:

         /**
          * @brief Real space scalar values
          */
         RTPScalar mRTP;

         /**
          * @brief Real space gradient field values
          */
         RTPField mGrad;
   };

   template<typename TSimType> inline PhysicalScalarFieldBase<TSimType>::PhysicalScalarFieldBase(SmartTruncation pTrunc, typename PhysicalScalarFieldBase<TSimType>::TransformType &transform)
      : PhysicalBase<TSimType>(pTrunc, transform), mRTP(pTrunc), mGrad(pTrunc)
   {
   }

   template<typename TSimType> inline const RTPScalar&  PhysicalScalarFieldBase<TSimType>::rtp() const
   {
      return this->mRTP;
   }

   template<typename TSimType> inline RTPScalar&  PhysicalScalarFieldBase<TSimType>::rRTP()
   {
      return this->mRTP;
   }

   template<typename TSimType> inline const RTPField&  PhysicalScalarFieldBase<TSimType>::grad() const
   {
      return this->mGrad;
   }

   template<typename TSimType> inline RTPField&  PhysicalScalarFieldBase<TSimType>::rGrad()
   {
      return this->mGrad;
   }

}

#endif // PHYSICALSCALARFIELDBASE_HPP
