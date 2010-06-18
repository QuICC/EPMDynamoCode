/** \file AssocLegendreBasis.cpp
 *  \brief Source of the implementation of a basis of Associated Legendre polynomials
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/AssocLegendreBasis.hpp"

// Project includes
//

namespace EPMDynamo {

   AssocLegendreBasis::AssocLegendreBasis(const ArrayI& parameters, const int nX, const int nL)
      :PolynomialBasis<AssociatedLegendreOperator>(parameters, nX, nL)
   {
      // Initialise the vector of basis polynomials
      this->initBasis();
   }

   AssocLegendreBasis::AssocLegendreBasis(const SmartTruncation pTrunc)
      :PolynomialBasis<AssociatedLegendreOperator>(pTrunc->local()->fdsh()->mArray(), pTrunc->sim()->hoz()->nTh(), pTrunc->sim()->hoz()->nL())
   {
      // Initialise the vector of basis polynomials
      this->initBasis();
   }

   void AssocLegendreBasis::initBasis()
   {
      // The actual size of the polynomial depends on its fixed parameter (nL = maxNL - m)
      int size = 0;
      size = this->polyN() - this->parameters()(0);
      // Full initiliasation of the first polynomial (including grid and weights computation)
      this->mBasis.push_back(AssociatedLegendreOperator(this->parameters()(0), this->gridN(), size));
      
      //Get pointer to grid points array
      SmartArray grid = at(0).pGrid();
      //Get pointer to weights array
      SmartArray weights = at(0).pWeights();

      // Initialise the remaining polynomials using the grid and weights pointers
      // In this way the grid and weights are only stored once. The other polynomials
      // simply share the pointer
      for(int i=1; i < this->basisN(); ++i)
      {
         // The actual size of the polynomial depends on its fixed parameter (nL = maxNL - m)
         size = this->polyN() - this->parameters()(i);
         this->mBasis.push_back(AssociatedLegendreOperator(this->parameters()(i), grid, size, weights));
      }
   }

}
