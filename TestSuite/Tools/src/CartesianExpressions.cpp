/** \file CartesianExpressions.cpp
 */

// System includes
//

// Class includes
//
#include "CartesianExpressions.hpp"

// Project includes
//
#include "General/MathConstants.hpp"

namespace EPMDynamo {

   CartesianExpressions::CartesianExpressions(SmartTruncation pTrunc)
      : mpTrunc(pTrunc), mNph(pTrunc->sim()->hoz()->nPh()),
         mRad(pTrunc->local()->rtp()->radGrid()),
         mSth(pTrunc->sim()->hoz()->sinTheta()),
         mCth(pTrunc->sim()->hoz()->cosTheta()),
         mPhi(pTrunc->sim()->hoz()->phGrid())
   {
   }

   Matrix CartesianExpressions::xMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);
      int th0 = this->trunc()->local()->rtp()->th0(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = mSth(th0 + th)*std::cos(mPhi(ph));
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::yMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);
      int th0 = this->trunc()->local()->rtp()->th0(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = mSth(th0 + th)*std::sin(mPhi(ph));
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::zMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);
      int th0 = this->trunc()->local()->rtp()->th0(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = mCth(th0 + th);
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::xToRMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);
      int th0 = this->trunc()->local()->rtp()->th0(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = mSth(th0 + th)*std::cos(mPhi(ph));
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::yToRMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);
      int th0 = this->trunc()->local()->rtp()->th0(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = mSth(th0 + th)*std::sin(mPhi(ph));
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::zToRMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);
      int th0 = this->trunc()->local()->rtp()->th0(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = mCth(th0 + th);
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::xToThMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);
      int th0 = this->trunc()->local()->rtp()->th0(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = mCth(th0 + th)*std::cos(mPhi(ph));
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::yToThMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);
      int th0 = this->trunc()->local()->rtp()->th0(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = mCth(th0 + th)*std::sin(mPhi(ph));
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::zToThMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);
      int th0 = this->trunc()->local()->rtp()->th0(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = -mSth(th0 + th);
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::xToPhMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = -std::sin(mPhi(ph));
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::yToPhMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);

      Matrix   tmp(mNph, nTh);

      for(int th = 0; th < nTh; th++)
      {
         for(int ph = 0; ph < mNph; ph++)
         {
            tmp(ph, th) = std::cos(mPhi(ph));
         }
      }

      return tmp;
   }

   Matrix CartesianExpressions::zToPhMat(const int n)
   {
      int nTh = this->trunc()->local()->rtp()->nTh(n);

      Matrix   tmp(mNph, nTh);

      tmp.setConstant(0.0);

      return tmp;
   }
}
