/** \file ImposedTorPolCurl.cpp
 *  \brief Implementation of Toroidal/Poloidal expanded imposed field
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalFields/Imposed/ImposedTorPolCurl.hpp"

// Project includes
//

namespace EPMDynamo {

   ImposedTorPolCurl::ImposedTorPolCurl(SmartTruncation pTrunc, ImposedTorPolCurl::TransformType &transform)
      : PhysicalTorPolCurlBase(pTrunc, transform), mTotalField(pTrunc), mImposedField(pTrunc)
   {
   }

   void ImposedTorPolCurl::initialiseZeros()
   {
      PhysicalTorPolCurlBase::initialiseZeros();

      this->mImposedField.initialiseZeros();
   }

}
