/** \file PhysicalRTPField.hpp
 *  \brief Base of the implementation of the RTP components of any physical field
 */

#ifndef PHYSICALRTPFIELD_HPP
#define PHYSICALRTPFIELD_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalBase.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralFields/RTPField.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the RTP components of any physical field
    *
    * This class contains the different element and representation that are required
    * for a generic physical field.
    */
   class PhysicalRTPField : public PhysicalBase
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

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
      protected:

      private:

         /**
          * @brief Real space field values
          */
         RTPField mRTP;
   };

   inline const RTPField&  PhysicalRTPField::rtp() const
   {
      return this->mRTP;
   }

   inline RTPField&  PhysicalRTPField::rRTP()
   {
      return this->mRTP;
   }
}

#endif // PHYSICALRTPFIELD_HPP
