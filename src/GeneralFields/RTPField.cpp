/** \file RTPField.cpp
 *  \brief Source of the RTP space vector field implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "GeneralFields/RTPField.hpp"

// Project includes
//

namespace EPMDynamo {

   RTPField::RTPField(const SmartTruncation pTrunc)
      : RTPBase(pTrunc), mRComponent(mpTrunc), mThComponent(mpTrunc), mPhComponent(mpTrunc)
   {
   }

   void RTPField::initialiseZeros()
   {
      // initialise the R component to zero
      this->rR().initialiseZeros();

      // initialise the Theta component to zero
      this->rTheta().initialiseZeros();

      // initialise the Phi component to zero
      this->rPhi().initialiseZeros();
   }

}
