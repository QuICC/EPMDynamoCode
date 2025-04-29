/** \file FDSHTruncation.cpp
 *  \brief Source of the implementation of the FDSH space truncation information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Core/FDSHTruncation.hpp"

// Project includes
//

namespace EPMDynamo {

   FDSHTruncation::FDSHTruncation(const ArrayI &r0, const ArrayI &nR, const ArrayI &mArray)
      : mR0(r0), mNr(nR), mNm(mArray.size()), mMArray(mArray)
   {
   }

   ArrayI   FDSHTruncation::nLArray(const int nL) const
   {
      ArrayI sizes(this->nM());

      for(int k=0; k < this->nM(); ++k)
      {
         sizes(k) = nL - this->mMArray(k);
      }

      return sizes;
   }

}
