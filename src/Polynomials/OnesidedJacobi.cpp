/** \file OnesidedJacobi.cpp
 *  \brief Source of the implementation of the onesided Jacobi polynomial
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/OnesidedJacobi.hpp"

// Project includes
//
#include "General/MathConstants.hpp"
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   OnesidedJacobi::OnesidedJacobi(const int pL, const int nX, const int nN)
      : WeightedPolynomial(pL, nX, nN, msMaxDiff), mAlpha(SimulationConstants::RadialGridTraits::ALPHA), mBeta(pL + SimulationConstants::RadialGridTraits::BETA), mBPoly(nN)
   {
      // Initialise everything
      this->initFull();
   }

   OnesidedJacobi::OnesidedJacobi(const int pL, const SmartArray grid, const int nN, const SmartArray weights)
      : WeightedPolynomial(pL, grid, nN, weights, msMaxDiff), mAlpha(SimulationConstants::RadialGridTraits::ALPHA), mBeta(pL + SimulationConstants::RadialGridTraits::BETA), mBPoly(nN)
   {
      // Initialise using the provided grid and weights
      this->initPartial();
   }

   void OnesidedJacobi::initFull()
   {
      // Compute the grid and weights
      this->computeWeightedGrid();

      // Compute the polynomials on the grid points
      this->computePolynomials();

      // Compute the polynomials on the grid boundary
      this->computeBoundaryPolynomials();

      // Compute the Weighted polynomials
      this->computeWPolynomials();

      // Compute the derivatives
      this->computeDerivatives();

      // Compute the boundary derivatives polynomials
      this->computeBoundaryDerivatives();

      // Compute the weighted derivatives
      this->computeWDerivatives();
   }

   void OnesidedJacobi::initPartial()
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

   void OnesidedJacobi::computeWeightedGrid()
   {
      // Initialise the grid generator
      SimulationConstants::RadialGridTraits::RadialGrid  generator(this->gridN(), GridGeneratorBase::RADIAL_GRID, SimulationConstants::isCSCSGrid());

      // get computed grid
      this->rGrid() = generator.grid();

      // get computed weights
      this->rWeights() = generator.weights();
   }

   void OnesidedJacobi::computePolynomials()
   {
      Array x(gridN());
      x.setConstant(-1.0);
      x.array() += 2.0*this->grid().array().pow(2);

      this->rPoly().row(0).setConstant(1.0);
      this->rPoly().row(1) = (this->recurrenceB(1)*x.transpose()).array() + this->recurrenceA(1);

      for(int n = 2; n < this->polyN(); ++n)
      {
         this->rPoly().row(n) = ((this->recurrenceB(n)*x.transpose()).array() + this->recurrenceA(n));
         this->rPoly().row(n).array() *= this->poly().row(n-1).array();
         this->rPoly().row(n) -= (this->recurrenceC(n)*this->poly().row(n-2));
      }
   }

   void OnesidedJacobi::computeBoundaryPolynomials()
   {
      this->rBPoly()(0) = 1.0;
      this->rBPoly()(1) = this->recurrenceB(1) + this->recurrenceA(1);

      for(int n = 2; n < polyN(); ++n)
      {
         this->rBPoly()(n) = this->recurrenceB(n) + this->recurrenceA(n);
         this->rBPoly()(n) *= this->bpoly()(n-1);
         this->rBPoly()(n) -= this->recurrenceC(n)*this->bpoly()(n-2);
      }
   }

   void OnesidedJacobi::computeWPolynomials()
   {
      for(int n = 0; n < this->polyN(); ++n)
      {
         this->rWPoly().col(n) = this->poly().row(n).transpose().array() * this->weights().array();
      }
   }

   void OnesidedJacobi::computeDerivatives()
   {
      // The derivatives should be computed from the highest to the lowest order
      // In this way the computation can be used several times
      //
      
      //Compute the second derivative
      this->computeSecondDerivative();

      //Compute the first derivative
      this->computeFirstDerivative();
   }

   void OnesidedJacobi::computeBoundaryDerivatives()
   {
      // Initiliase diff vector
      this->initBDiff();

      // The derivatives should be computed from the highest to the lowest order
      // In this way the computation can be used several times
      
      //Compute the second derivative
      this->computeBoundarySecondDerivative();

      //Compute the first derivative
      this->computeBoundaryFirstDerivative();
   }

   void OnesidedJacobi::computeDerivativeBase(Matrix& rMat, const int k) const
   {
      Array x(gridN());
      x.setConstant(-1.0);
      x.array() += 2.0*grid().array().pow(2);

      EPMFloat dalpha = this->alpha() + static_cast<EPMFloat>(k);
      EPMFloat dbeta = this->beta() + static_cast<EPMFloat>(k);

      for(int n=0; n < k; ++n)
      {
         rMat.row(n).setConstant(0.0);
      }
      rMat.row(k).setConstant(1.0);
      rMat.row(k+1) = (this->recurrenceB(1, dalpha, dbeta)*x.transpose()).array() + this->recurrenceA(1,dalpha,dbeta);

      for(int n = k+2; n < this->polyN(); ++n)
      {
         rMat.row(n) = ((this->recurrenceB(n-k, dalpha, dbeta)*x.transpose()).array() + this->recurrenceA(n-k,dalpha,dbeta));
         rMat.row(n).array() *= rMat.row(n-1).array();
         rMat.row(n) -= (this->recurrenceC(n-k,dalpha,dbeta)*rMat.row(n-2));
      }

      for(int n = k; n < this->polyN(); ++n)
      {
         for(int m=1; m <= k; ++m)
         {
            rMat.row(n) *= (static_cast<EPMFloat>(n)+this->alpha()+this->beta()+static_cast<EPMFloat>(m));
         }
      }
   }

   void OnesidedJacobi::computeBoundaryDerivativeBase(Array& rMat, const int k) const
   {
      EPMFloat dalpha = this->alpha() + static_cast<EPMFloat>(k);
      EPMFloat dbeta = this->beta() + static_cast<EPMFloat>(k);

      for(int n=0; n < k; ++n)
      {
         rMat(n) = 0.0;
      }
      rMat(k) = 1.0;
      rMat(k+1) = this->recurrenceB(1, dalpha, dbeta) + this->recurrenceA(1,dalpha,dbeta);

      for(int n = k+2; n < polyN(); ++n)
      {
         rMat(n) = this->recurrenceB(n-k, dalpha, dbeta) + this->recurrenceA(n-k,dalpha,dbeta);
         rMat(n) *= rMat(n-1);
         rMat(n) -= this->recurrenceC(n-k,dalpha,dbeta)*rMat(n-2);
      }

      for(int n = k; n < polyN(); ++n)
      {
         for(int m=1; m <= k; ++m)
         {
            rMat(n) *= (static_cast<EPMFloat>(n)+this->alpha()+this->beta()+static_cast<EPMFloat>(m));
         }
      }
   }

   void OnesidedJacobi::initBDiff()
   {
      for(int i=0; i < this->msMaxDiff; ++i)
      {
         this->mBDiff.push_back(Array(this->polyN()));
      }
   }

   void OnesidedJacobi::computeFirstDerivative()
   {
      // This supposes that the second derivative already intialised it
      for(int n = 1; n < this->polyN(); ++n)
      {
         this->rDiff(1).row(n).array() *= 2.0*this->grid().array();
      }
   }

   void OnesidedJacobi::computeBoundaryFirstDerivative()
   {
      // This supposes that the second derivative already intialised it
      for(int n = 1; n < this->polyN(); ++n)
      {
         this->rBDiff(1)(n) *= 2.0;
      }
   }

   void OnesidedJacobi::computeSecondDerivative()
   {
      // Initialise with base values of first derivatives
      this->computeDerivativeBase(this->rDiff(1), 1);

      // Initialise with base values of first derivatives
      this->computeDerivativeBase(this->rDiff(2), 2);

      for(int n = 1; n < polyN(); ++n)
      {
         this->rDiff(2).row(n).array() *= 4.0*this->grid().array().pow(2).transpose();
         this->rDiff(2).row(n) += 2.0*this->diff(1).row(n);
      }
   }

   void OnesidedJacobi::computeBoundarySecondDerivative()
   {
      // Initialise with base values of first derivatives
      this->computeBoundaryDerivativeBase(this->rBDiff(1), 1);

      // Initialise with base values of first derivatives
      this->computeBoundaryDerivativeBase(this->rBDiff(2), 2);

      for(int n = 1; n < this->polyN(); ++n)
      {
         this->rBDiff(2)(n) *= 4.0;
         this->rBDiff(2)(n) += 2.0*this->bdiff(1)(n);
      }
   }

   void OnesidedJacobi::computeWDerivatives()
   {
      for(int i=1; i <= this->msMaxDiff; ++i)
      {
         for(int n = 0; n < this->polyN(); ++n)
         {
            this->rWDiff(i).col(n).array() = this->diff(i).row(n).transpose().array() * this->weights().array();
         }
      }
   }

   EPMFloat OnesidedJacobi::recurrenceA(const int n) const
   {
      return this->recurrenceA(n, this->alpha(), this->beta());
   }

   EPMFloat OnesidedJacobi::recurrenceB(const int n) const
   {
      return this->recurrenceB(n, this->alpha(), this->beta());
   }

   EPMFloat OnesidedJacobi::recurrenceC(const int n) const
   {
      return this->recurrenceC(n, this->alpha(), this->beta());
   }

   EPMFloat OnesidedJacobi::recurrenceA(const int n, const EPMFloat alpha, const EPMFloat beta) const
   {
      if(n > 1)
      {
         EPMFloat factorA;
         EPMFloat albe = alpha + beta;
         EPMFloat a2_b2 = alpha*alpha - beta*beta;
         EPMFloat rn = static_cast<EPMFloat>(n);

         factorA = (2*rn + albe - 1)*a2_b2/(2*rn*(rn + albe)*(2*rn + albe - 2.0));

         return factorA;
      } else if(n == 1)
      {
         return (alpha-beta)/2.0;
      } else
      {
         return 0.0;
      }
   }

   EPMFloat OnesidedJacobi::recurrenceB(const int n, const EPMFloat alpha, const EPMFloat beta) const
   {
      EPMFloat albe = alpha + beta;

      if(n > 1)
      {
         EPMFloat factorB;
         EPMFloat rn = static_cast<EPMFloat>(n);

         factorB = (2*rn + albe - 1)*(2*rn + albe)/(2*rn*(rn + albe));

         return factorB;
      } else if (n == 1)
      {
         return albe/2.0 + 1.0;
      } else
      {
         return 0.0;
      }

   }

   EPMFloat OnesidedJacobi::recurrenceC(const int n, const EPMFloat alpha, const EPMFloat beta) const
   {
      if(n > 1)
      {
         EPMFloat factorC;
         EPMFloat albe = alpha + beta;
         EPMFloat rn = static_cast<EPMFloat>(n);

         factorC = (2*rn + albe)*(rn + alpha - 1.0)*(rn + beta - 1.0)/(rn*(rn + albe)*(2*rn + albe - 2.0));

         return factorC;
      } else
      {
         return 0.0;
      }
   }
}
