/** \file WeightedPolynomial.cpp
 *  \brief Source of the implementation of a general weighted polynomial
 */

// System includes
//

// Class include
//
#include "Polynomials/WeightedPolynomial.hpp"

// Project includes
//

namespace EPMDynamo {

   WeightedPolynomial::WeightedPolynomial(const int param, const int nX, const int nN, const int maxDiff)
      : Polynomial(param, nX, nN, maxDiff), mpWeights(new Array(nX)), mWPoly(nX, nN)
   {
      // Initialise the weighted derivatives
      this->initWDiff();
   }

   WeightedPolynomial::WeightedPolynomial(const int param, const SmartArray grid, const int nN, const SmartArray weights, const int maxDiff)
      : Polynomial(param, grid, nN, maxDiff), mpWeights(weights), mWPoly(grid->size(), nN)
   {
      // Initialise the weighted derivatives
      this->initWDiff();
   }

   void WeightedPolynomial::initWDiff()
   {
      // Initialise the vector of matrix for the derivatives
      for(int i=0; i < maxDiff(); ++i)
      {
         this->mWDiff.push_back(Matrix(this->gridN(), this->polyN()));
      }
   }
}
