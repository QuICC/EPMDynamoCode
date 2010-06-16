/** \file LegendreGrid.cpp
 *  \brief Source of the Legendre grid generator
 */

// System includes
//

// External includes
//

// Class include
//
#include "GridGenerators/LegendreGrid.hpp"

// Project includes
//
#include "Polynomials/QuadratureRule.hpp"
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   LegendreGrid::LegendreGrid(const int gridN, const int type, bool isCSCS)
      : GridGeneratorBase(gridN, type, isCSCS)
   {
      // Generate grid and weights
      this->initialise();
   }

   void LegendreGrid::initialise()
   {
      // Generate a grid and weights
      this->computeWeightedGrid();

      if(this->mType == GridGeneratorBase::RADIAL_GRID)
      {
         // Convert to radial grid
         this->convertToRadial();

         // Sort the grid in increasing ordering
         this->sortGrid();
      } 
   }

   void LegendreGrid::addCSCS(int &left, int &right)
   {
      if(this->mType == GridGeneratorBase::X_GRID)
      {
         this->addCSCSTheta(left, right);
      } else if(this->mType == GridGeneratorBase::RADIAL_GRID)
      {
         this->addCSCSRadial(left, right);
      }
   }

   void LegendreGrid::computeWeightedGrid(int left, int right)
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
      QuadratureRule    quadRule(diag, subdiag, 2.0);

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

   void LegendreGrid::computeQuadratureMatrix(Array& diag, Array& subdiag)
   {
      int ptsN = diag.size();

      // The diagonal elements are zero
      diag.setZero();

      // Loop over the subdiagonal terms
      for(int n=1; n < ptsN; ++n)
      {
         subdiag(n-1) = std::sqrt(static_cast<DynamoFloat>(n*n)/static_cast<DynamoFloat>(4*n*n-1));
      }
   }

   void LegendreGrid::addCSCSRadial(int &left, int &right)
   {
      // modify grid shift factors
      left++;

      // Set additional r = 1 point
      this->mGrid(this->gridN() - 1) = 1.0;
      this->mWeights(this->gridN() - 1) = 0.0;
   }

   void LegendreGrid::addCSCSTheta(int &left, int &right)
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
