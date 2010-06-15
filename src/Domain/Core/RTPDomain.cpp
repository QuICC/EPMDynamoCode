/** \file RTPDomain.cpp
 *  \brief Source of the implementation of the RTP space domain information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Core/RTPDomain.hpp"

// Project includes
//

namespace EPMDynamo {

   RTPDomain::RTPDomain(const int r0, const int nR, const ArrayI &th0, const ArrayI &nTh)
      : RTPTruncation(r0, nR, th0, nTh)
   {
   }

   void RTPDomain::setRadialGrid(SmartArray rad)
   {
      if(rad->size() == nR())
      {
         this->mpRadGrid = rad;
      } else
      {
         this->mpRadGrid = SmartArray(new Array(this->nR()));
         for(int n=0; n < this->nR(); ++n)
         {
            (*this->mpRadGrid)(n) = (*rad)(this->r0() + n);
         }
      }
   }

   void RTPDomain::setThetaGrid(SmartArray cTh, SmartArray sTh)
   {
      this->mpCosTh = cTh;
      this->mpSinTh = sTh;
   }

}
