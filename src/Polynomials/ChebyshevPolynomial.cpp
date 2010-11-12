/** \file ChebyshevPolynomial.cpp
 *  \brief Source of the implementation of the Chebyshev polynomial
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/ChebyshevPolynomial.hpp"

// Project includes
//
#include "General/MathConstants.hpp"
#include "Simulations/SimulationConstants.hpp"
#include "GridGenerators/ChebyshevGrid.hpp"

namespace EPMDynamo {

   ChebyshevPolynomial::ChebyshevPolynomial(const int pL, const int nX, const int nN)
      : WeightedPolynomial(pL, nX, nN, msMaxDiff), mBPoly(nN)
   {
      // Initialise everything
      this->initFull();
   }

   ChebyshevPolynomial::ChebyshevPolynomial(const int pL, const SmartArray grid, const int nN, const SmartArray weights)
      : WeightedPolynomial(pL, grid, nN, weights, msMaxDiff), mBPoly(nN)
   {
      // Initialise using the provided grid and weights
      this->initPartial();
   }

   void ChebyshevPolynomial::initFull()
   {
      // Compute the grid points
      this->computeWeightedGrid();

      // Compute the polynomials on the grid points
      this->computePolynomials();

      // Compute the polynomials on the grid boundary
      this->computeBoundaryPolynomials();

      // Compute the derivatives on the grid points
      this->computeWPolynomials();

      // Compute the weighted polynomials
      this->computeDerivatives();

      // Compute the weighted polynomials
      this->computeBoundaryDerivatives();

      // Compute the weighted derivatives
      this->computeWDerivatives();
   }

   void ChebyshevPolynomial::initPartial()
   {
      // Compute the polynomials on the grid points
      this->computePolynomials();

      // Compute the polynomials on the boundary
      this->computeBoundaryPolynomials();

      // Compute the weighted polynomials
      this->computeWPolynomials();

      // Compute the derivatives on the grid points
      this->computeDerivatives();

      // Compute the derivatives on the grid points
      this->computeBoundaryDerivatives();

      // Compute the weighted derivatives
      this->computeWDerivatives();
   }

   void ChebyshevPolynomial::computeWeightedGrid()
   {
      // Initialise the grid generator
      ChebyshevGrid  generator(this->gridN(), GridGeneratorBase::RADIAL_GRID, SimulationConstants::isCSCSGrid());

      // get computed grid
      this->rGrid() = generator.grid();

      // get computed weights
      this->rWeights() = generator.weights();
   }

   void ChebyshevPolynomial::computePolynomials()
   {
      Array x(this->gridN());
      x.setConstant(-1.0);
      x.array() += 2.0*this->grid().array().pow(2);

      for(int i = 0; i < x.size(); ++i)
      {
         x(i) = std::acos(x(i));
      }

      this->rPoly().row(0).setConstant(1.0);

      for(int n = 1; n < this->polyN(); ++n)
      {
         this->rPoly().row(n).array() = (static_cast<EPMFloat>(n)*x).transpose().array().cos();
      }
   }

   void ChebyshevPolynomial::computeBoundaryPolynomials()
   {
      this->rBPoly()(0) = 1.0;
      this->rBPoly()(1) = this->recurrenceB(1) + this->recurrenceA(1);

      for(int n = 2; n < this->polyN(); ++n)
      {
         this->rBPoly()(n) = this->recurrenceB(n) + this->recurrenceA(n);
         this->rBPoly()(n) *= this->bpoly()(n-1);
         this->rBPoly()(n) -= this->recurrenceC(n)*this->bpoly()(n-2);
      }
   }

   void ChebyshevPolynomial::computeWPolynomials()
   {
      for(int n = 0; n < this->polyN(); ++n)
      {
         this->rWPoly().col(n) = (this->poly().row(n).transpose().array() * this->weights().array() * this->normaliseW(n)).matrix();
      }
   }

   void ChebyshevPolynomial::computeDerivatives()
   {
      //Compute the second derivative
      this->computeSecondDerivative();

      //Compute the first derivative
      this->computeFirstDerivative();
   }

   void ChebyshevPolynomial::computeBoundaryDerivatives()
   {
      // Initiliase diff vector
      this->initBDiff();
      
      //Compute the second derivative
      this->computeBoundarySecondDerivative();

      //Compute the first derivative
      this->computeBoundaryFirstDerivative();
   }

   void ChebyshevPolynomial::initBDiff()
   {
      for(int i=0; i < this->msMaxDiff; ++i)
      {
         this->mBDiff.push_back(Array(this->polyN()));
      }
   }

   void ChebyshevPolynomial::computeFirstDerivative()
   {
      Array x(this->gridN());
      x.setConstant(-1.0);
      x.array() += 2.0*this->grid().array().pow(2);

      for(int i = 0; i < x.size(); ++i)
      {
         x(i) = std::acos(x(i));
      }

      this->rDiff(1).row(0).setConstant(0.0);
      this->rDiff(1).row(1).setConstant(1.0);

      EPMFloat rn;
      for(int n = 2; n < this->polyN(); ++n)
      {
         rn = static_cast<EPMFloat>(n);

         // Compute dr part
         this->rDiff(1).row(n).array() = 4.0*(0.5*x).array().cos().transpose();

         // Compute the trigonometric derivative part
         this->rDiff(1).row(n).array() *= (rn*((rn*x).array().sin()/x.array().sin())).transpose();
      }
   }

   void ChebyshevPolynomial::computeBoundaryFirstDerivative()
   {
      this->rBDiff(1)(0) = 0.0;

      EPMFloat rn;
      for(int n = 1; n < this->polyN(); ++n)
      {
         rn = static_cast<EPMFloat>(n);

         this->rBDiff(1)(n) = 4.0*(rn*rn);
      }
   }

   void ChebyshevPolynomial::computeSecondDerivative()
   {
      Array x(this->gridN());
      x.setConstant(-1.0);
      x.array() += 2.0*this->grid().array().pow(2);

      for(int i = 0; i < x.size(); ++i)
      {
         x(i) = std::acos(x(i));
      }

      this->rDiff(2).row(0).setConstant(0.0);
      this->rDiff(2).row(1).setConstant(4.0);

      EPMFloat rn;
      for(int n = 2; n < this->polyN(); ++n)
      {
         rn = static_cast<EPMFloat>(n);

         this->rDiff(2).row(n).array() = (-2.0*rn*(rn*x).array().cos()).transpose();

         this->rDiff(2).row(n).array() += (((0.5*x).array().sin()/(0.5*x).array().cos())*(rn*x).array().sin()).transpose();

         this->rDiff(2).row(n).array() *= (4.0*rn*(((0.5*x).array().sin()/(0.5*x).array().cos())/x.array().sin())).transpose();
      }
   }

   void ChebyshevPolynomial::computeBoundarySecondDerivative()
   {
      this->rBDiff(2)(0) = 0.0;

      EPMFloat rn;
      for(int n = 1; n < this->polyN(); ++n)
      {
         rn = static_cast<EPMFloat>(n);

         this->rBDiff(2)(n) = 4.0/3.0*rn*rn*(4.0*rn*rn-1.0);
      }
   }

   void ChebyshevPolynomial::computeWDerivatives()
   {
      for(int i=1; i <= this->msMaxDiff; ++i)
      {
         for(int n = 0; n < this->polyN(); ++n)
         {
            this->rWDiff(i).col(n).array() = this->diff(i).row(n).transpose().array() * this->weights().array()* this->normaliseW(n);
         }
      }
   }

   EPMFloat ChebyshevPolynomial::recurrenceA(const int n) const
   {
     return 0.0;
   }

   EPMFloat ChebyshevPolynomial::recurrenceB(const int n) const
   {
      return 2.0;
   }

   EPMFloat ChebyshevPolynomial::recurrenceC(const int n) const
   {
      return 1.0;
   }

   EPMFloat ChebyshevPolynomial::normaliseW(const int n)
   {
      return 1.0/this->weights().dot(this->poly().row(n).array().pow(2).matrix());
   }
}
