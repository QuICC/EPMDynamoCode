/** \file PhysicalRTPScalar.hpp
 *  \brief Base of the implementation of the RTP component of a physical scalar
 */

#ifndef PHYSICALRTPSCALAR_HPP
#define PHYSICALRTPSCALAR_HPP

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

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the RTP component of a physical scalar
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalRTPScalar : public PhysicalBase<TSimType>
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
         PhysicalRTPScalar(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalRTPScalar() {};

         /**
          * @brief Get the RTP scalar values
          */
         const RTPScalar&   rtp() const;

         /**
          * @brief Set the RTP scalar values
          */
         RTPScalar&   rRTP();
         
      protected:

      private:

         /**
          * @brief Real space scalar values
          */
         RTPScalar mRTP;
   };

   template<typename TSimType> inline PhysicalRTPScalar<TSimType>::PhysicalRTPScalar(SmartTruncation pTrunc, typename PhysicalRTPScalar<TSimType>::TransformType &transform)
      : PhysicalBase<TSimType>(pTrunc, transform), mRTP(pTrunc)
   {
   }

   template<typename TSimType> inline const RTPScalar&  PhysicalRTPScalar<TSimType>::rtp() const
   {
      return this->mRTP;
   }

   template<typename TSimType> inline RTPScalar&  PhysicalRTPScalar<TSimType>::rRTP()
   {
      return this->mRTP;
   }

}

#endif // PHYSICALRTPSCALAR_HPP
