/** \file PhysicalRTPScalar.hpp
 *  \brief Base of the implementation of the RTP component of a physical scalar
 */

#ifndef PHYSICALRTPSCALAR_HPP
#define PHYSICALRTPSCALAR_HPP

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
#include "GeneralScalars/RTPScalar.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the RTP component of a physical scalar
    */
   class PhysicalRTPScalar : public PhysicalBase
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

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

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
      protected:

      private:

         /**
          * @brief Real space scalar values
          */
         RTPScalar mRTP;
   };

   inline const RTPScalar&  PhysicalRTPScalar::rtp() const
   {
      return this->mRTP;
   }

   inline RTPScalar&  PhysicalRTPScalar::rRTP()
   {
      return this->mRTP;
   }

}

#endif // PHYSICALRTPSCALAR_HPP
