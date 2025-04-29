/** \file PhysicalQSTCurlBase.cpp
 *  \brief Implementation of the base of the QST expanded physical field
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalFields/Fields/PhysicalQSTCurlBase.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalQSTCurlBase::PhysicalQSTCurlBase(SmartTruncation pTrunc, PhysicalQSTCurlBase::TransformType &transform)
      : PhysicalRTPFieldCurl(pTrunc, transform), mNeedCurlTransform(0), mPerturbation(pTrunc)
   {
   }

   void PhysicalQSTCurlBase::initialiseZeros()
   {
      PhysicalRTPFieldCurl::initialiseZeros();

      this->mPerturbation.initialiseZeros();
   }

}
