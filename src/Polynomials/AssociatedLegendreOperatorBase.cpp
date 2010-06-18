/** \file AssociatedLegendreOperatorBase.cpp
 *  \brief Source of the base of the implementation of an operator set based on the associated Legendre polynomial
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/AssociatedLegendreOperatorBase.hpp"

// Project includes
//
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   AssociatedLegendreOperatorBase::AssociatedLegendreOperatorBase(const int pM, const int nX, const int nL)
      : AssociatedLegendrePolynomial(pM, nX, nL), mSin_1Theta(nX), mSll_1(nL), mSin_1Ones(), mSin_1Poly()
   {
      // Compute the corresponding 1/sin(theta) values
      this->computeSin_1Theta();

      // Compute the corresponding 1/sqrt(l(l+1)) values
      this->computeSll_1();

      // Treat the CSCS case differently (has 2 extraneous grid points)
      if(SimulationConstants::isCSCSGrid())
      {
         // Setup the ones and P/Sin for CSCS grids
         this->computeSpecialOperators();
      }
   }

   AssociatedLegendreOperatorBase::AssociatedLegendreOperatorBase(const int pM, const SmartArray grid, const int nL, const SmartArray weights)
      : AssociatedLegendrePolynomial(pM, grid, nL, weights), mSin_1Theta(grid->size()), mSll_1(nL), mSin_1Ones(), mSin_1Poly()
   {
      // Compute the corresponding 1/sin(theta) values
      this->computeSin_1Theta();

      // Compute the corresponding 1/sqrt(l(l+1)) values
      this->computeSll_1();

      // Treat the CSCS case differently (has 2 extraneous grid points)
      if(SimulationConstants::isCSCSGrid())
      {
         // Setup the ones and P/Sin for CSCS grids
         this->computeSpecialOperators();
      }
   }

   void AssociatedLegendreOperatorBase::computeSin_1Theta()
   {
      // Treat the CSCS case differently (has 2 extraneous grid points)
      if(SimulationConstants::isCSCSGrid())
      {
         // Compute the array of 1.0/sin(theta) coordinates from grid points theta values
         this->mSin_1Theta(0) = 0.0;
         for(int i=1; i < this->gridN() - 1; ++i)
         {
            this->mSin_1Theta(i) = 1.0/this->rSinTheta()(i);
         }
         this->mSin_1Theta(this->gridN()-1) = 0.0;
      } else
      {
         // Compute the array of 1.0/sin(theta) coordinates from grid points theta values
         for(int i=0; i < this->gridN(); ++i)
         {
            this->mSin_1Theta(i) = 1.0/this->rSinTheta()(i);
         }
      }
   }

   void AssociatedLegendreOperatorBase::computeSll_1()
   {
      // Compute the array of 1.0/sqrt(l(l+1))
      if(this->m() == 0)
      {
         this->mSll_1(0) = 0.0;
      } else {
         this->mSll_1(0) = 1.0/std::sqrt(static_cast<EPMFloat>(this->m()*(this->m()+1)));
      }

      for(int i= 1, l = this->m() + 1 ; i < this->polyN(); ++i, ++l)
      {
         this->mSll_1(i) = 1.0/std::sqrt(static_cast<EPMFloat>(l*(l+1)));
      }
   }

   void AssociatedLegendreOperatorBase::computeSpecialOperators()
   {
      // Set the special array for 1/Sin(theta)
      this->mSin_1Ones = this->mSin_1Theta;
      this->mSin_1Ones(0) = 1.0;
      this->mSin_1Ones(this->mSin_1Ones.size()-1) = 1.0;

      // Set the special Poly/Sin matrix
      this->mSin_1Poly = this->poly();

      // Set the field values to the phi=0 coordinate system value
      for(int i = 0; i < this->diff(1).rows(); ++i)
      {
//         this->mSin_1Poly(i,0) = -this->diff(1)(0,0);
//         this->mSin_1Poly(i,this->mSin_1Ones.size()-1) = this->diff(1)(0,this->mSin_1Ones.size()-1);
      }
   }
}
