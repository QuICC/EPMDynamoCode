/** \file AssociatedLegendrePolynomial.cpp
 *  \brief Source of the implementation of the Associated Legendre polynomial
 *  This file contains the definitions of the AssociatedLegendrePolynomial class methods
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/AssociatedLegendrePolynomial.hpp"

// Project includes
//
#include "General/EPMException.hpp"
#include "Polynomials/QuadratureRule.hpp"
#include "Simulations/SimulationConstants.hpp"
#include "GridGenerators/LegendreGrid.hpp"

namespace EPMDynamo {

   AssociatedLegendrePolynomial::AssociatedLegendrePolynomial(const int pM, const int nX, const int nL)
      : WeightedPolynomial(pM, nX, nL, this->msMaxDiff), mpTheta(new Array(nX)), mpSinTheta(new Array(nX))
   {
      // Initialise everything
      this->initFull();
   }

   AssociatedLegendrePolynomial::AssociatedLegendrePolynomial(const int pM, const SmartArray grid, const int nL, const SmartArray weights)
      : WeightedPolynomial(pM, grid, nL, weights, msMaxDiff), mpTheta(new Array(grid->size())), mpSinTheta(new Array(grid->size()))
   {
      // Initialise using the provided grid and weights
      this->initPartial();
   }

   void AssociatedLegendrePolynomial::initFull()
   {
      // Compute the weights and grid points
      this->computeWeightedGrid();

      // Compute the grid points
      this->computeThetaGrids();

      // Compute the polynomials on the grid points
      this->computePolynomials();

      // Compute the derivatives on the grid points
      this->computeDerivatives();

      // Compute the weighted polynomials
      this->computeWPolynomials();

      // Compute the weighted derivatives
      this->computeWDerivatives();
   }

   void AssociatedLegendrePolynomial::initPartial()
   {
      // Compute the grid points
      this->computeThetaGrids();

      // Compute the polynomials on the grid points
      this->computePolynomials();

      // Compute the derivatives on the grid points
      this->computeDerivatives();

      // Compute the weighted polynomials
      this->computeWPolynomials();

      // Compute the weighted derivatives
      this->computeWDerivatives();
   }

   void AssociatedLegendrePolynomial::computeWeightedGrid()
   {
      // Initialise the grid generator
      LegendreGrid  generator(this->gridN(), GridGeneratorBase::X_GRID, SimulationConstants::isCSCSGrid());

      // get computed grid
      this->rGrid() = generator.grid();

      // get computed weights
      this->rWeights() = generator.weights();
   }

   void AssociatedLegendrePolynomial::computeThetaGrids()
   {
      // Compute the corresponding theta values
      this->computeTheta();

      // Compute the corresponding sin(theta) values
      this->computeSinTheta();
   }

   void AssociatedLegendrePolynomial::computeTheta()
   {
      // Compute the array of theta coordinates from grid points (cos(theta))
      for(int i=0; i < this->gridN(); ++i)
      {
         this->rTheta()(i) = std::acos(this->grid()(i));
      }
   }

   void AssociatedLegendrePolynomial::computeSinTheta()
   {
      // Compute the array of sin(theta) coordinates from grid points theta values
      for(int i=0; i < this->gridN(); ++i)
      {
         this->rSinTheta()(i) = std::sin(this->theta()(i));
      }
   }

   EPMFloat AssociatedLegendrePolynomial::normaliseW(const int l) const
   {
      // Set basis value of the normalisation factor (due to schmidt normalisation)
      EPMFloat norm = (2.0*static_cast<EPMFloat>(l) + 1.0)/4.0;

      // An additional factor of 2 is required for m=0
      if(m() == 0)
      {
         return 2.0*norm;
      } else
      {
         return norm;
      }
   }

   void AssociatedLegendrePolynomial::computePolynomials()
   {
      // Compute the Plm polynomials with fixed m
      this->computePlmBasis();
   }

   void AssociatedLegendrePolynomial::computeWPolynomials()
   {
      // Initialise the weighted polynomials including normalisation factor
      for(int n = 0; n < this->polyN(); ++n)
      {
         this->rWPoly().col(n) = this->poly().row(n).transpose().array() * this->weights().array() * this->normaliseW(n+m());
      }
   }

   void AssociatedLegendrePolynomial::computeDerivatives()
   {
      // Compute the polynomials first derivatives with fixed m
      this->computeDPlmBasis();
   }

   void AssociatedLegendrePolynomial::computeWDerivatives()
   {
      // Initialise the weighted polynomials derivatives including normalisation factor
      for(int i=1; i <= this->msMaxDiff; ++i)
      {
         for(int n = 0; n < this->polyN(); ++n)
         {
            this->rWDiff(i).col(n) = this->diff(i).row(n).transpose().array() * this->weights().array() * this->normaliseW(n+m());
         }
      }
   }

   EPMFloat AssociatedLegendrePolynomial::recurrenceA(const int n) const
   {
      return 0.0;
   }

   EPMFloat AssociatedLegendrePolynomial::recurrenceB(const int n) const
   {
      if(n > 1)
      {
         EPMFloat factorB;
         EPMFloat rn = static_cast<EPMFloat>(n);

         factorB = (2*rn - 1)/rn;

         return factorB;
      } else if (n == 1)
      {
         return 1.0;
      } else
      {
         return 0.0;
      }

   }

   EPMFloat AssociatedLegendrePolynomial::recurrenceC(const int n) const
   {
      if(n > 1)
      {
         EPMFloat factorC;
         EPMFloat rn = static_cast<EPMFloat>(n);

         factorC = (rn - 1.0)/rn;

         return factorC;
      } else
      {
         return 0.0;
      }
   }

   void AssociatedLegendrePolynomial::computePmm(const int m, const int idx)
   {
      if (m < 0)
      {
         throw EPMException("AssociatedLegendrePolynomial::computePmm", "m < 0");
      } else if (m == 0)
      {
         this->rPoly().row(idx).setConstant(1.0);
      } else
      {
         EPMFloat r2i;
         EPMFloat factor = std::sqrt(2.0);

         for (int i=1; i <= m; i++)
         {
            r2i = static_cast<EPMFloat>(2*i);
            factor *= -std::sqrt((r2i - 1.0)/r2i);
         }
         this->rPoly().row(idx) = this->sinTheta().transpose().array().pow(m);
         this->rPoly().row(idx) *= factor;
      }
   }

   void AssociatedLegendrePolynomial::computePmm1(const int m, const int idx)
   {
      if (m < 0)
      {
         throw EPMException("AssociatedLegendrePolynomial::computePmm1", "m < 0");
      }
      this->rPoly().row(idx) = this->grid().transpose();
      this->rPoly().row(idx) *= std::sqrt(2.0*static_cast<EPMFloat>(m) + 1.0);
      this->rPoly().row(idx).array() *= this->poly().row(idx-1).array();
   }

   void AssociatedLegendrePolynomial::computePlm(int l, const int m, int idx)
   {
      if (l < 2)
      {
         throw EPMException("AssociatedLegendrePolynomial::computePlm", "l < 2");
      }

      if ((m < 0) || (m > l))
      {
         throw EPMException("AssociatedLegendrePolynomial::computePlm", "m < 0 or m > l");
      }

      if (l == m)
      {
         throw EPMException("AssociatedLegendrePolynomial::computePlm", "l = m");
      }

      EPMFloat rm = static_cast<EPMFloat>(m);
      EPMFloat rl = static_cast<EPMFloat>(l);

      int maxL = l + this->polyN() - 2;

      for (; l < maxL; ++l, ++idx)
      {
         this->rPoly().row(idx) = (2.0*rl - 1.0)*this->grid();
         this->rPoly().row(idx).array() *= this->poly().row(idx-1).array();

         this->rPoly().row(idx) -= this->poly().row(idx-2)*std::sqrt((rl + rm - 1.0)*(rl - rm - 1.0));

         this->rPoly().row(idx) /= std::sqrt((rl + rm)*(rl - rm));

         rl += 1.0;
      }
   }

   void AssociatedLegendrePolynomial::computeDPmm(const int m, const int idx)
   {
      if (m < 0) 
      {
         throw EPMException("AssociatedLegendrePolynomial::computeDPmm", "m < 0");
      } else if (m == 0)
      {
         this->rDiff(1).row(idx).setZero();
      } else
      {

         EPMFloat r2i;
         EPMFloat factor = std::sqrt(2.0)*static_cast<EPMFloat>(m);
         for (int i=1; i <= m; i++)
         {
            r2i = static_cast<EPMFloat>(2*i);
            factor *= -std::sqrt((r2i-1.0)/r2i);
         }

         this->rDiff(1).row(idx) = this->sinTheta().transpose().array().pow(m-1);
         this->rDiff(1).row(idx).array() *= this->grid().transpose().array();
         this->rDiff(1).row(idx) *= factor;
      
      }
   }

   void AssociatedLegendrePolynomial::computeDPmm1(const int m, const int idx)
   {

      if (m < 0) 
      {
         throw EPMException("AssociatedLegendrePolynomial::computeDPmm1", "m < 0");
      }

      EPMFloat rm = static_cast<EPMFloat>(m);

      this->rDiff(1).row(idx) = (this->diff(1).row(idx-1).array()*this->grid().transpose().array() - this->poly().row(idx-1).array()*this->sinTheta().transpose().array());
      this->rDiff(1).row(idx) *= std::sqrt(2.0*rm + 1.0);

   }

   void AssociatedLegendrePolynomial::computeDPlm(int l, const int m, int idx)
   {
      if (l < 2)
      {
         throw EPMException("AssociatedLegendrePolynomial::computeDPlm", "l < 2");
      }

      if ((m < 0) || (m > l))
      {
         throw EPMException("AssociatedLegendrePolynomial::computeDPlm", "m < 0 or m > l");
      }

      if (l == m)
      {
         throw EPMException("AssociatedLegendrePolynomial::computeDPlm", "l = 0");
      }

      EPMFloat rm = static_cast<EPMFloat>(m);
      EPMFloat rl = static_cast<EPMFloat>(l);

      int maxL = l + polyN() - 2;

      for (; l < maxL; ++l, ++idx)
      {
         this->rDiff(1).row(idx) = this->diff(1).row(idx-1).array()*this->grid().transpose().array();
         this->rDiff(1).row(idx) -= (this->poly().row(idx-1).array()*this->sinTheta().transpose().array()).matrix();
         this->rDiff(1).row(idx) *= (2.0*rl - 1.0);

         this->rDiff(1).row(idx) -= this->diff(1).row(idx-2)*std::sqrt((rl + rm - 1.0)*(rl - rm - 1.0));

         this->rDiff(1).row(idx) /= std::sqrt((rl + rm)*(rl - rm));

         rl += 1.0;
      }
   }

   void AssociatedLegendrePolynomial::computePlmBasis()
   {
      // Compute P_m^m
      //
      this->computePmm(m());

      if(polyN() > 1)
      {
         // Compute P_{m+1}^m
         //
         this->computePmm1(m());

         if(this->polyN() > 2)
         {
            // Compute P_l^m
            //
            this->computePlm(m()+2, m());
         }
      }
   }

   void AssociatedLegendrePolynomial::computeDPlmBasis()
   {
      // Compute dP_m^m
      //
      this->computeDPmm(m());

      if(this->polyN() > 1)
      {
         // Compute dP_{m+1}^m
         //
         this->computeDPmm1(m());

         if(this->polyN() > 2)
         {
            // Compute dP_l^m
            //
            this->computeDPlm(m()+2, m());
         }
      }
   }
}
