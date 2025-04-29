/** \file PhysicalRTPScalarGradient.cpp
 *  \brief Base of the implementation of the RTP component of a physical scalar with gradient
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalRTPScalarGradient.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalRTPScalarGradient::PhysicalRTPScalarGradient(SmartTruncation pTrunc, PhysicalRTPScalarGradient::TransformType &transform)
      : PhysicalRTPScalar(pTrunc, transform), mGrad(pTrunc)
   {
   }

   void  PhysicalRTPScalarGradient::initialiseZeros()
   {
      PhysicalRTPScalar::initialiseZeros();

      this->mGrad.initialiseZeros();
   }

}
