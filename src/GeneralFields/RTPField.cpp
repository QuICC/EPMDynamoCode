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

   RTPField::RTPField(const SmartSTrunc pSTrunc)
      : RTPBase(pSTrunc), mRComponent(mpSTrunc), mThComponent(mpSTrunc), mPhComponent(mpSTrunc)
   {
   }

}
