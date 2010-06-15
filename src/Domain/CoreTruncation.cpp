/** \file CoreTruncation.cpp
 *  \brief Source of the core truncation information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/CoreTruncation.hpp"

// Project includes
//

namespace EPMDynamo {

   CoreTruncation::CoreTruncation(const int rR0, const int rNr, const ArrayI &rTh0, const ArrayI &rNth, const ArrayI &fR0, const ArrayI &fNr, const ArrayI &fM, const ArrayI &sL, const std::vector<ArrayI> &sMs)
      : CoreTruncationBase(fR0, fNr, fM, sL, sMs), mpRTP(new RTPTruncation(rR0, rNr, rTh0, rNth))
   {
   }
}
