/** \file PhysicalRTPField.hpp
 *  \brief Base of the implementation of the RTP components of any physical field
 */

#ifndef PHYSICALRTPFIELD_HPP
#define PHYSICALRTPFIELD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalBase.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralFields/RTPField.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the RTP components of any physical field
    *
    * This class contains the different element and representation that are required
    * for a generic physical field.
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalRTPField : public PhysicalBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /**
         * @brief Constructs the field elements and representations for the given
         * spectral truncation
         *
         * @param pTrunc Smart pointer to truncation information
         * @param transform Transform object
         */
         PhysicalRTPField(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Simple empty destructor
         */
         virtual ~PhysicalRTPField() {};

         /**
          * @brief Get the RTP field values
          */
         const RTPField&   rtp() const;

         /**
          * @brief Set the RTP field values
          */
         RTPField&   rRTP();
         
      protected:

      private:

         /**
          * @brief Real space field values
          */
         RTPField mRTP;
   };

   template<typename TSimType> inline PhysicalRTPField<TSimType>::PhysicalRTPField(SmartTruncation pTrunc, typename PhysicalRTPField<TSimType>::TransformType &transform)
      : PhysicalBase<TSimType>(pTrunc, transform), mRTP(pTrunc)
   {
   }

   template<typename TSimType> inline const RTPField&  PhysicalRTPField<TSimType>::rtp() const
   {
      return this->mRTP;
   }

   template<typename TSimType> inline RTPField&  PhysicalRTPField<TSimType>::rRTP()
   {
      return this->mRTP;
   }
}

#endif // PHYSICALRTPFIELD_HPP
