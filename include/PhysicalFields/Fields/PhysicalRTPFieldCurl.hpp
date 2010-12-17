/** \file PhysicalRTPFieldCurl.hpp
 *  \brief Base of the implementation of the RTP components of any physical field
 *  with curl
 */

#ifndef PHYSICALRTPFIELDCURL_HPP
#define PHYSICALRTPFIELDCURL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalRTPField.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the RTP components of any physical field
    * with curl
    *
    * This class contains the different element and representation that are required
    * for a generic physical field.
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalRTPFieldCurl : public PhysicalRTPField<TSimType>
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
         PhysicalRTPFieldCurl(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Simple empty destructor
         */
         virtual ~PhysicalRTPFieldCurl() {};

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
          * @brief Real space curl of the field
          */
         RTPField mCurl;
   };

   template<typename TSimType> inline PhysicalRTPFieldCurl<TSimType>::PhysicalRTPFieldCurl(SmartTruncation pTrunc, typename PhysicalRTPFieldCurl<TSimType>::TransformType &transform)
      : PhysicalRTPField<TSimType>(pTrunc, transform), mCurl(pTrunc)
   {
   }

   template<typename TSimType> inline const RTPField&  PhysicalRTPFieldCurl<TSimType>::curl() const
   {
      return this->mCurl;
   }

   template<typename TSimType> inline RTPField&  PhysicalRTPFieldCurl<TSimType>::rCurl()
   {
      return this->mCurl;
   }

}

#endif // PHYSICALRTPFIELDCURL_HPP
