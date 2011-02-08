/** \file RadialDomain.cpp
 *  \brief Source of the implementation of the radial domain information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Simulation/RadialDomain.hpp"

// Project includes
//
#include "GridGenerators/JacobiGrid.hpp"
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   RadialDomain::RadialDomain(const int maxN, const int nR)
      : RadialTruncation(maxN, nR), mMaxL(0), mpRSll(new Array(nR))
   {
   }

   EPMFloat RadialDomain::alpha() const
   {
      return -0.5;
   }

   EPMFloat RadialDomain::beta(int l) const
   {
      return static_cast<EPMFloat>(l)-0.5;
   }

   void RadialDomain::setRadSll()
   {
      // Create the grid generator
      int gridN = static_cast<int>(std::ceil(static_cast<EPMFloat>(3*this->maxN())/2.0));
      JacobiGrid   generator(gridN, GridGeneratorBase::RADIAL_GRID, SimulationConstants::isCSCSGrid());

      // Create storage for the jacobiGrid
      Array jacobiGrid(this->nN());

      // Loop over all harmonic degrees starting from the highest
      for(int l = this->mMaxL; l > 0; --l)
      {
         // Recompute the grid for given parameters
         generator.recompute(this->alpha(), this->beta(l));

         // Get the computed grid
         jacobiGrid = generator.grid();

         // Loop over the radiuses
         int r = 0;
         while(jacobiGrid(0) > this->radGrid()(r))
         {
            r++;
            (*this->mpRSll)(r) = std::sqrt(static_cast<EPMFloat>(l*(l+1)));
         }

         // Special case for l=0
         (*this->mpRSll)(0) = std::sqrt(2.0);

      }
   }

}
