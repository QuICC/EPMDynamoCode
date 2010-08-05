/** \file WorlandPolynomial.cpp
 *  \brief Source of the implementation of the Worland polynomial
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/WorlandPolynomial.hpp"

// Project includes
//
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   WorlandPolynomial::WorlandPolynomial(const int pL, const int nX, const int nN)
      : OnesidedJacobi(pL, nX, nN)
   {
      // Initialise everything
      this->initFull();
   }

   WorlandPolynomial::WorlandPolynomial(const int pL, const SmartArray grid, const int nN, const SmartArray weights)
      : OnesidedJacobi(pL, grid, nN, weights)
   {
      // Initialise using the provided grid and weights
      this->initPartial();
   }

   void WorlandPolynomial::initFull()
   {
      // Order is important. Derivatives should be computed first from highest order to 
      // lowest ...
      
      // Compute the derivatives
      this->computeDerivatives();

      // Compute the boundary derivatives polynomials
      this->computeBoundaryDerivatives();

      // Compute the weighted derivatives
      this->computeWDerivatives();

      // ... Then only compute the polynomials. 
     
      // Compute the polynomials on the grid points
      this->computePolynomials();

      // Compute the Weighted polynomials
      this->computeWPolynomials();

      // Correct the polynomials for special cases (for example CSCS output)
      this->correctPolynomials();
   }

   void WorlandPolynomial::initPartial()
   {
      // Order is important. Derivatives should be computed first from highest order to 
      // lowest ...
      
      // Compute the derivatives
      this->computeDerivatives();

      // Compute the boundary derivatives polynomials
      this->computeBoundaryDerivatives();

      // Compute the weighted derivatives
      this->computeWDerivatives();

      // ... Then only compute the polynomials. 
      
      // Compute the polynomials on the grid points
      this->computePolynomials();

      // Compute the Weighted polynomials
      this->computeWPolynomials();

      // Correct the polynomials for special cases (for example CSCS output)
      this->correctPolynomials();
   }

   EPMFloat WorlandPolynomial::normaliseW(const int n)
   {
      return 1.0/this->weights().dot(this->poly().row(n).cwise().pow(2));
   }

   void WorlandPolynomial::computePolynomials()
   {
      if(this->l() != 0)
      {
         Array   rfactor = this->grid().cwise().pow(l()); 
         for(int n = 0; n < this->polyN(); ++n)
         {
            this->rPoly().row(n).cwise() *= rfactor.transpose();
         }
      }
   }

   void WorlandPolynomial::computeDerivatives()
   {
      // The derivatives should be computed from the highest to the lowest order
      // In this way the computation can be used several times
      
      //Compute the second derivative
      this->computeSecondDerivative();

      //Compute the first derivative
      this->computeFirstDerivative();
   }

   void WorlandPolynomial::computeBoundaryDerivatives()
   {
      // The derivatives should be computed from the highest to the lowest order
      // In this way the computation can be used several times
       
      //Compute the second derivative
      this->computeBoundarySecondDerivative();

      //Compute the first derivative
      this->computeBoundaryFirstDerivative();
   }

   void WorlandPolynomial::computeFirstDerivative()
   {
      Array rfactor = static_cast<EPMFloat>(l())*this->grid().cwise().pow(l()-1);
      Array drfactor = this->grid().cwise().pow(l());

      for(int n = 0; n < this->polyN(); ++n)
      {
         this->rDiff(1).row(n) = drfactor.transpose().cwise()*this->diff(1).row(n);
         this->rDiff(1).row(n) += rfactor.transpose().cwise()*this->poly().row(n);
      }
   
   }

   void WorlandPolynomial::computeBoundaryFirstDerivative()
   {
      EPMFloat rfactor = static_cast<EPMFloat>(l());

      for(int n = 0; n < polyN(); ++n)
      {
         this->rBDiff(1)(n) = this->bdiff(1)(n);
         this->rBDiff(1)(n) += rfactor*this->bpoly()(n);
      }
   
   }

   void WorlandPolynomial::computeSecondDerivative()
   {
      Array rfactor = static_cast<EPMFloat>(l()*(l()-1))*this->grid().cwise().pow(l()-2);
      Array drfactor = static_cast<EPMFloat>(2*l())*this->grid().cwise().pow(l()-1);
      Array ddrfactor = this->grid().cwise().pow(l());

      for(int n = 0; n < this->polyN(); ++n)
      {
         this->rDiff(2).row(n) = ddrfactor.transpose().cwise()*this->diff(2).row(n);
         this->rDiff(2).row(n) += drfactor.transpose().cwise()*this->diff(1).row(n);
         this->rDiff(2).row(n) += rfactor.transpose().cwise()*this->poly().row(n);
      }
   
   }

   void WorlandPolynomial::computeBoundarySecondDerivative()
   {
      EPMFloat rfactor = static_cast<EPMFloat>(this->l()*(this->l()-1));
      EPMFloat drfactor = static_cast<EPMFloat>(2*this->l());

      for(int n = 0; n < this->polyN(); ++n)
      {
         this->rBDiff(2)(n) = this->bdiff(2)(n);
         this->rBDiff(2)(n) += drfactor*this->bdiff(1)(n);
         this->rBDiff(2)(n) += rfactor*this->bpoly()(n);
      }
   
   }

   void WorlandPolynomial::computeWPolynomials()
   {
      for(int n = 0; n < this->polyN(); ++n)
      {
         this->rWPoly().col(n) = this->poly().row(n).transpose().cwise() * this->weights() * this->normaliseW(n);
      }
   }

   void WorlandPolynomial::computeWDerivatives()
   {
      for(int i=1; i <= this->msMaxDiff; ++i)
      {
         for(int n = 0; n < this->polyN(); ++n)
         {
            this->rWDiff(i).col(n) = this->diff(i).row(n).transpose().cwise() * this->weights() * this->normaliseW(n);
         }
      }
   }

   void WorlandPolynomial::correctPolynomials()
   {
      // Treat the CSCS case differently (has 2 extraneous grid points)
      if(SimulationConstants::isCSCSGrid())
      {
         // Set special values for r=0 and l < 1
         if(this->l() < 1)
         {
            this->rDiff(1).col(0).setConstant(0.0);
            this->rWDiff(1).row(0).setConstant(0.0);
         }

         // Set special values for r=0 and l < 2
         if(this->l() < 2)
         {
            this->rDiff(2).col(0).setConstant(0.0);
            this->rWDiff(2).row(0).setConstant(0.0);
         }
      }
   }
}
