/** \file PhysicalTorPolCurlBase.cpp
 *  \brief Base of the implementation of Toroidal/Poloidal expanded field field
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
#include "PhysicalFields/Fields/PhysicalTorPolCurlBase.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalTorPolCurlBase::PhysicalTorPolCurlBase(SmartTruncation pTrunc, PhysicalTorPolCurlBase::TransformType &transform)
      : PhysicalRTPFieldCurl(pTrunc, transform), mNeedCurlTransform(0), mPerturbation(pTrunc)
   {
   }

   void  PhysicalTorPolCurlBase::initialiseZeros()
   {
      PhysicalRTPFieldCurl::initialiseZeros();

      this->mPerturbation.initialiseZeros();
   }

}
