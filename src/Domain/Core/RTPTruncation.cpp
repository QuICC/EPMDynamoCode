/** \file RTPTruncation.cpp
 *  \brief Source of the implementation of the RTP space truncation information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Core/RTPTruncation.hpp"

// Project includes
//

namespace EPMDynamo {

   RTPTruncation::RTPTruncation(const int r0, const int nR, const ArrayI &th0, const ArrayI &nTh)
      : mR0(r0), mNr(nR), mTh0(th0), mNth(nTh)
   {
   }

}
