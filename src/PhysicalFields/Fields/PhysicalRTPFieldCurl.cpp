/** \file PhysicalRTPFieldCurl.cpp
 *  \brief Base of the implementation of the RTP components of any physical field
 *  with curl
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalFields/Fields/PhysicalRTPFieldCurl.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalRTPFieldCurl::PhysicalRTPFieldCurl(SmartTruncation pTrunc, PhysicalRTPFieldCurl::TransformType &transform)
      : PhysicalRTPField(pTrunc, transform), mCurl(pTrunc)
   {
   }

   void  PhysicalRTPFieldCurl::initialiseZeros()
   {
      PhysicalRTPField::initialiseZeros();

      this->mCurl.initialiseZeros();
   }

}
