/** \file Polynomial.cpp
 *  \brief Source of the implementation of a general polynomial
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/Polynomial.hpp"

// Project includes
//

namespace EPMDynamo {

   Polynomial::Polynomial(const int param, const int nX, const int nN, const int maxDiff)
      : mParam(param), mGridN(nX), mPolyN(nN), mMaxDiff(maxDiff), mpGrid(new Array(nX)), mPoly(nN, nX)
   {
      // Initialise the storage for the derivatives
      this->initDiff();
   }

   Polynomial::Polynomial(const int param, const SmartArray grid, const int nN, const int maxDiff)
      : mParam(param), mGridN(grid->size()), mPolyN(nN), mMaxDiff(maxDiff), mpGrid(grid), mPoly(nN, mGridN)
   {
      // Initialise the storage for the derivatives
      this->initDiff();
   }

   void Polynomial::initDiff()
   {
      // Initialise the vector of Matrices for the derivatives
      for(int i=0; i < this->mMaxDiff; ++i)
      {
         this->mDiff.push_back(Matrix(this->mPolyN, this->mGridN));
      }
   }
}
