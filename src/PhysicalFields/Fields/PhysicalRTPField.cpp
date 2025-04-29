/** \file PhysicalRTPField.cpp
 *  \brief Base of the implementation of the RTP components of any physical field
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalRTPField.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalRTPField::PhysicalRTPField(SmartTruncation pTrunc, PhysicalRTPField::TransformType &transform)
      : PhysicalBase(pTrunc, transform), mRTP(pTrunc)
   {
   }

   void PhysicalRTPField::initialiseZeros()
   {
      this->mRTP.initialiseZeros();
   }
}
