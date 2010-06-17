/** \file FFTransformBase.cpp
 *  \brief Source of the base of the implementation of the FFT transform
 */

// System includes
//

// External includes
//

// Class include
//
#include "Transforms/FFT/FFTransformBase.hpp"

// Project includes
//

namespace EPMDynamo {

   FFTransformBase::FFTransformBase(SmartTruncation pTrunc)
      : mSize(pTrunc->sim()->hoz()->nPh()), mZRow(mSize/2 + 1), mHowmany(0), mReSize(0), mZSize(0), mScale(1.0)
   {
      // Get size of transforms
      this->mHowmany = (pTrunc->local()->rtp()->nThArray().sum());

      // Init size array
      this->initSizes();

      // Initialiation of the scaling factor (FFTW does unnormalised transforms)
      this->mScale = 1.0/static_cast<double>(this->mSize);
   }

   void FFTransformBase::initSizes()
   {
      // Set the size of the real array
      this->mReSize = this->mSize*this->mHowmany;

      // Set the size of the complex array
      this->mZSize = this->mZRow*this->mHowmany;
   }

}
