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

}
