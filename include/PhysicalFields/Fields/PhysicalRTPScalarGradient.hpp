/** \file PhysicalRTPScalarGradient.hpp
 *  \brief Base of the implementation of the RTP component of a physical scalar with gradient
 */

#ifndef PHYSICALRTPSCALARGRADIENT_HPP
#define PHYSICALRTPSCALARGRADIENT_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalRTPScalar.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralFields/RTPField.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the RTP component of a physical scalar with gradient
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalRTPScalarGradient : public PhysicalRTPScalar<TSimType>
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
         PhysicalRTPScalarGradient(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalRTPScalarGradient() {};

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
          * @brief Real space gradient field values
          */
         RTPField mGrad;
   };

   template<typename TSimType> inline PhysicalRTPScalarGradient<TSimType>::PhysicalRTPScalarGradient(SmartTruncation pTrunc, typename PhysicalRTPScalarGradient<TSimType>::TransformType &transform)
      : PhysicalRTPScalar<TSimType>(pTrunc, transform), mGrad(pTrunc)
   {
   }

   template<typename TSimType> inline const RTPField&  PhysicalRTPScalarGradient<TSimType>::grad() const
   {
      return this->mGrad;
   }

   template<typename TSimType> inline RTPField&  PhysicalRTPScalarGradient<TSimType>::rGrad()
   {
      return this->mGrad;
   }

}

#endif // PHYSICALRTPSCALARGRADIENT_HPP
