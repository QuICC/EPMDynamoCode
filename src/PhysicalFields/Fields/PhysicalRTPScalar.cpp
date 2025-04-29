/** \file PhysicalRTPScalar.cpp
 *  \brief Base of the implementation of the RTP component of a physical scalar
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalFields/Fields/PhysicalRTPScalar.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalRTPScalar::PhysicalRTPScalar(SmartTruncation pTrunc, PhysicalRTPScalar::TransformType &transform)
      : PhysicalBase(pTrunc, transform), mRTP(pTrunc)
   {
   }

   void  PhysicalRTPScalar::initialiseZeros()
   {
      this->mRTP.initialiseZeros();
   }

}
