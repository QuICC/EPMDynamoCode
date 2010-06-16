/** \file ChebyshevGrid.cpp
 *  \brief Source of the Chebyshev grid generator
 */

// System includes
//

// External includes
//

// Class include
//
#include "GridGenerators/ChebyshevGrid.hpp"

// Project includes
//
#include "General/MathConstants.hpp"
#include "General/EPMException.hpp"
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   ChebyshevGrid::ChebyshevGrid(const int gridN, const int type, bool isCSCS)
      : GridGeneratorBase(gridN, type, isCSCS)
   {
      // Generate grid and weights
      this->initialise();
   }

   void ChebyshevGrid::initialise()
   {
      // Compute grid and weights
      this->computeWeightedGrid();

      if(this->mType == GridGeneratorBase::RADIAL_GRID)
      {
         // Convert to radial grid
         this->convertToRadial();

         // Sort the grid in increasing ordering
         this->sortGrid();
      } else
      {
         throw EPMException("ChebyshevGrid::initialise", "Unknown grid type!");
      }
   }

   void ChebyshevGrid::computeWeightedGrid(int left, int right)
   {
      // Handle CSCS output special grid
      if(this->mIsCSCS)
      {
         // modify shifting due to CSCS points
         this->addCSCSRadial(left, right);
      }

      // Compute (partial) grid
      this->computeGrid(left, right);

      // Compute (partial) grid weights
      this->computeWeights(left, right);
   }

   void ChebyshevGrid::computeGrid(const int left, const int right)
   {
      int xN = this->gridN() - left - right;
      DynamoFloat f_tmp;

      for(int i=0; i < xN; ++i)
      {
         f_tmp = (MathConstants::PI/2.0)*static_cast<DynamoFloat>((2*(i+1)-1))/static_cast<DynamoFloat>(xN);

         this->mGrid(i + left) = std::cos(f_tmp);
      }
   }

   void ChebyshevGrid::computeWeights(const int left, const int right)
   {
      int xN = this->gridN() - left - right;

      for(int i=0; i < xN; ++i)
      {
         this->mWeights(i + left)  = MathConstants::PI/static_cast<DynamoFloat>(2*xN);
      }
   }

   void ChebyshevGrid::addCSCSRadial(int &left, int &right)
   {
      // modify grid shift factors
      left++;
      right++;

      // Set additional r = 0 point
      this->mGrid(0) = -1.0;
      this->mWeights(0) = 0.0;

      // Set additional r = 1 point
      this->mGrid(this->gridN() - 1) = 1.0;
      this->mWeights(this->gridN() - 1) = 0.0;
   }

}
