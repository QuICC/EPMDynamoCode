/** \file RTPScalar.cpp
 *  \brief Source of the RTP scalar field data type
 */

// System includes
//

// External includes
//

// Class include
//
#include "GeneralScalars/RTPScalar.hpp"

// Project includes
//

namespace EPMDynamo {

   RTPScalar::RTPScalar(const SmartTruncation pTrunc)
      : RTPBase(pTrunc)
   {
      this->initShells();
   }

   void RTPScalar::initShells()
   {
      int nR = this->nR();
      for(int n=0; n < nR; ++n)
      {
         this->mValues.push_back(SphericalShell(this->nPh(), this->nTh(n)));
      }
   }

   void RTPScalar::initialiseZeros()
   {
      int nR = this->nR();
      for(int n=0; n < nR; ++n)
      {
         this->rShell(n).setConstant(0.0);
      }
   }

}
