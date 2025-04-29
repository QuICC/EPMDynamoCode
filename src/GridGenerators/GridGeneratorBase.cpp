/** \file GridGeneratorBase.cpp
 *  \brief Source of the base for the grid generators
 */

// System includes
//
#include <map>

// External includes
//

// Class include
//
#include "GridGenerators/GridGeneratorBase.hpp"

// Project includes
//
#include "General/MathConstants.hpp"
#include "General/EPMException.hpp"

namespace EPMDynamo {

   GridGeneratorBase::GridGeneratorBase(const int gridN, const int type, bool isCSCS)
      : mGrid(gridN), mWeights(gridN), mType(type), mIsCSCS(isCSCS)
   {
   }

   void GridGeneratorBase::convertToRadial()
   {
      // compute x + 1
      this->mGrid.array() += 1.0;

      // compute (x+1)/2
      this->mGrid /= 2.0;

      // compute ((x+1)/2)^(1/2)
      this->mGrid = this->mGrid.array().sqrt();
   }

   void GridGeneratorBase::sortGrid()
   {
      // Create a map to sort elements
      std::map<EPMFloat, EPMFloat> sorter;

      // fill map with grid points /weights pairs
      for(int i = 0; i < this->gridN(); ++i)
      {
         sorter.insert(std::make_pair(this->mGrid(i), this->mWeights(i)));
      }

      // Check that no point got lost
      if(sorter.size() != static_cast<unsigned int>(this->gridN()))
      {
//         throw EPMException("GridGeneratorBase::sortGrid", "Lost grid points during conversion!");
      }

      std::map<EPMFloat, EPMFloat>::const_iterator  it;

      // Replace grid point values with reorder version
      int i = 0;
      for(it = sorter.begin(); it != sorter.end(); ++it, ++i)
      {
         this->mGrid(i) = it->first;
         this->mWeights(i) = it->second;
      }
   }

}
