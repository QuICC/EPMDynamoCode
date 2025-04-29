/** \file CoreTruncationBase.cpp
 *  \brief Source of the base of the core/remote truncation information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/CoreTruncationBase.hpp"

// Project includes
//

namespace EPMDynamo {

   CoreTruncationBase::CoreTruncationBase(const ArrayI &fR0, const ArrayI &fNr, const ArrayI &fM, const ArrayI &sL, const std::vector<ArrayI> &sMs)
      : mpFDSH(new FDSHTruncation(fR0, fNr, fM)), mpSpec(new SpecTruncation(sL, sMs))
   {
   }

}
