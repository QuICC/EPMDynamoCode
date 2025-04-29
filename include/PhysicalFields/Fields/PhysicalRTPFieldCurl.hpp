/** \file PhysicalRTPFieldCurl.hpp
 *  \brief Base of the implementation of the RTP components of any physical field
 *  with curl
 */

#ifndef PHYSICALRTPFIELDCURL_HPP
#define PHYSICALRTPFIELDCURL_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalRTPField.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the RTP components of any physical field
    * with curl
    *
    * This class contains the different element and representation that are required
    * for a generic physical field.
    */
   class PhysicalRTPFieldCurl : public PhysicalRTPField
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

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

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
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

   inline const RTPField&  PhysicalRTPFieldCurl::curl() const
   {
      return this->mCurl;
   }

   inline RTPField&  PhysicalRTPFieldCurl::rCurl()
   {
      return this->mCurl;
   }

}

#endif // PHYSICALRTPFIELDCURL_HPP
