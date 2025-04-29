/** \file PhysicalRTPScalarGradient.hpp
 *  \brief Base of the implementation of the RTP component of a physical scalar with gradient
 */

#ifndef PHYSICALRTPSCALARGRADIENT_HPP
#define PHYSICALRTPSCALARGRADIENT_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalRTPScalar.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralFields/RTPField.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the RTP component of a physical scalar with gradient
    */
   class PhysicalRTPScalarGradient : public PhysicalRTPScalar
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

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();
         
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

   inline const RTPField&  PhysicalRTPScalarGradient::grad() const
   {
      return this->mGrad;
   }

   inline RTPField&  PhysicalRTPScalarGradient::rGrad()
   {
      return this->mGrad;
   }

}

#endif // PHYSICALRTPSCALARGRADIENT_HPP
