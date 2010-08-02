/** \file PhysicalFieldBase.hpp
 *  \brief Base of the implementation of any physical field
 */

#ifndef PHYSICALFIELDBASE_HPP
#define PHYSICALFIELDBASE_HPP

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
    * \brief Base of the implementation of any physical field
    *
    * This class contains the different element and representation that are required
    * for a generic physical field.
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalFieldBase : public PhysicalBase<TSimType>
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
         PhysicalFieldBase(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Simple empty destructor
         */
         virtual ~PhysicalFieldBase() {};

         /**
          * @brief Get the RTP field values
          */
         const RTPField&   rtp() const;

         /**
          * @brief Set the RTP field values
          */
         RTPField&   rRTP();

         /**
          * @brief Get values of the curl of the RTP field
          */
         const RTPField&   curl() const;

         /**
          * @brief Compute RTP values of the curl of the field
          *
          * @param step Current step in a multistep transform
          */
         virtual void curlTransform(const int step) = 0;
         
      protected:

         /**
          * @brief Set values of the curl of the RTP field
          */
         RTPField&   rCurl();

      private:

         /**
          * @brief Real space field values
          */
         RTPField mRTP;

         /**
          * @brief Real space curl of the field
          */
         RTPField mCurl;
   };

   template<typename TSimType> inline PhysicalFieldBase<TSimType>::PhysicalFieldBase(SmartTruncation pTrunc, typename PhysicalFieldBase<TSimType>::TransformType &transform)
      : PhysicalBase<TSimType>(pTrunc, transform), mRTP(pTrunc), mCurl(pTrunc)
   {
   }

   template<typename TSimType> inline const RTPField&  PhysicalFieldBase<TSimType>::rtp() const
   {
      return this->mRTP;
   }

   template<typename TSimType> inline RTPField&  PhysicalFieldBase<TSimType>::rRTP()
   {
      return this->mRTP;
   }

   template<typename TSimType> inline const RTPField&  PhysicalFieldBase<TSimType>::curl() const
   {
      return this->mCurl;
   }

   template<typename TSimType> inline RTPField&  PhysicalFieldBase<TSimType>::rCurl()
   {
      return this->mCurl;
   }

}

#endif // PHYSICALFIELDBASE_HPP
