/** \file FDSHOrderScalarBase.cpp
 *  \brief Source of the base of the harmonic order ordered data type
 */

// System includes
//

// External includes
//

// Class include
//
#include "GeneralScalars/FDSHOrderScalarBase.hpp"

// Project includes
//

namespace EPMDynamo {

   FDSHOrderScalarBase::FDSHOrderScalarBase(const SmartTruncation pTrunc)
      : mpTrunc(pTrunc)
   {
      this->mLSizes = this->mpTrunc->local()->fdsh()->nLArray(this->mpTrunc->sim()->hoz()->nL());
   }

}
