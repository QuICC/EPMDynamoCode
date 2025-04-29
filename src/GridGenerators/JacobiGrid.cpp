/** \file JacobiGrid.cpp
 *  \brief Source of the Jacobi grid generator
 */

// System includes
//

// External includes
//

// Class include
//
#include "GridGenerators/JacobiGrid.hpp"

// Project includes
//
#include "General/EPMException.hpp"
#include "Polynomials/QuadratureRule.hpp"
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   JacobiGrid::JacobiGrid(const int gridN, const int type, bool isCSCS)
      : GridGeneratorBase(gridN, type, isCSCS), mAlpha(-0.5), mBeta(-0.5)
   {
      // Generate grid and weights
      this->initialise();
   }

   void JacobiGrid::recompute(EPMFloat alpha, EPMFloat beta)
   {
      // Set new alpha value
      this->mAlpha = alpha;

      // Set new beta value
      this->mBeta = beta;

      // Recompute the grid
      this->initialise();
   }

   void JacobiGrid::initialise()
   {
      // Generate a grid and weights
      this->computeWeightedGrid();

      if(this->mType == GridGeneratorBase::RADIAL_GRID)
      {
         // Convert to radial grid
         this->convertToRadial();

         // Sort the grid in increasing ordering
         this->sortGrid();
      } else
      {
         throw EPMException("JacobiGrid::initialise", "Unknown grid type!");
      }
   }

   void JacobiGrid::addCSCS(int &left, int &right)
   {
      if(this->mType == GridGeneratorBase::X_GRID)
      {
         this->addCSCSTheta(left, right);
      } else if(this->mType == GridGeneratorBase::RADIAL_GRID)
      {
         this->addCSCSRadial(left, right);
      }
   }

   void JacobiGrid::computeWeightedGrid(int left, int right)
   {
      // Handle CSCS output special grid
      if(this->mIsCSCS)
      {
         // modify shifting due to CSCS points
         this->addCSCS(left, right);
      }

      // Allow for only a partial grid computation using the left and right offsets
      int ptsN = this->gridN() - left - right;

      // Storage for the diagonal and subdiagonal values
      Array    diag(ptsN);
      Array    subdiag(ptsN - 1);

      // Compute the diagonal and subdiagonal elements of the matrix
      this->computeQuadratureMatrix(diag, subdiag);

      // Create a quadrature rule object
      QuadratureRule    quadRule(diag, subdiag, 1.0);

      // Check for a partial grid computation
      if(left + right == 0)
      {
         // Compute the grid points and the related weights
         quadRule.computeGridAndWeights(this->mGrid, this->mWeights);
      } else
      {
         // Create temporary storage
         Array tmpGrid(ptsN);
         Array tmpWeights(ptsN);

         // Compute the grid points and the related weights
         quadRule.computeGridAndWeights(tmpGrid, tmpWeights);

         // store the obtained partial grid at right place
         for(int i = 0; i < ptsN; ++i)
         {
            this->mGrid(i + left) = tmpGrid(i);
            this->mWeights(i + left) = tmpWeights(i);
         }
      }
   }

   void JacobiGrid::computeQuadratureMatrix(Array& diag, Array& subdiag)
   {
      int ptsN = diag.size();

      // Compute the diagonal elements
      EPMFloat dn;

      // Treat the n=1 case separately due to numerical issues
      diag(0) = (this->mBeta - this->mAlpha);
      if(diag(0) != 0.0)
      {
         diag(0) /= (2.0 + this->mAlpha + this->mBeta);
      }

      // Create the matrix entries for n > 1
      for(int n=2; n <= ptsN; ++n)
      {
         dn = static_cast<EPMFloat>(n);
         diag(n-1) = (this->mBeta*this->mBeta - this->mAlpha*this->mAlpha);
         if(diag(n-1) != 0.0)
         {
            diag(n-1) /= (2.0*dn + this->mAlpha + this->mBeta - 2.0)*(2.0*dn + this->mAlpha + this->mBeta);
         }
      }

      // Loop over the subdiagonal terms
      
      // Treat the n=1 case separately due to numerical issues
      subdiag(0) = (1.0 + this->mAlpha)*(1.0 + this->mBeta);
      if(subdiag(0) != 0.0)
      {
         subdiag(0) /= (2.0 + this->mAlpha + this->mBeta)*(2.0 + this->mAlpha + this->mBeta)*(3.0 + this->mAlpha + this->mBeta);
      }
      subdiag(0) = 2.0*std::sqrt(subdiag(0));

      // Create the matrix entries for n > 1
      for(int n=2; n < ptsN; ++n)
      {
         dn = static_cast<EPMFloat>(n);
         subdiag(n-1) = 4.0*dn*(dn + this->mAlpha)*(dn + this->mBeta)*(dn + this->mAlpha + this->mBeta);
         if(subdiag(n-1) != 0.0)
         {
            subdiag(n-1) /= (2.0*dn + this->mAlpha + this->mBeta)*(2.0*dn + this->mAlpha + this->mBeta)*(2.0*dn + this->mAlpha + this->mBeta + 1.0)*(2.0*dn + this->mAlpha + this->mBeta - 1.0);
         }
         subdiag(n-1) = std::sqrt(subdiag(n-1));
      }

      // Special case
      if(this->mAlpha == -0.5 && this->mBeta == -0.5)
      {
         subdiag(0) = 0.5;
      }
   }

   void JacobiGrid::addCSCSRadial(int &left, int &right)
   {
      // modify grid shift factors
      left++;
      right++;

      // Set additional r = 0 point
      this->mGrid(0) = 0.0;
      this->mWeights(0) = 0.0;

      // Set additional r = 1 point
      this->mGrid(this->gridN() - 1) = 1.0;
      this->mWeights(this->gridN() - 1) = 0.0;
   }

   void JacobiGrid::addCSCSTheta(int &left, int &right)
   {
      // modify grid shift factors
      left++;
      right++;

      // Set additional theta = pi point
      this->mGrid(0) = -1.0;
      this->mWeights(0) = 0.0;

      // Set additional theta = 0 point
      this->mGrid(this->gridN() - 1) = 1.0;
      this->mWeights(this->gridN() - 1) = 0.0;
   }

}
