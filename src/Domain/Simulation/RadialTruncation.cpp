/** \file RadialTruncation.cpp
 *  \brief Source of the implementation of the radial truncation information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Simulation/RadialTruncation.hpp"

// Project includes
//
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   RadialTruncation::RadialTruncation(const int maxN, const int nR)
      : mMaxN(maxN), mNn(maxN+1), mNr(nR)
   {
      // Correct grid sizes if required (for example with CSCS grid)
      this->correctGridSizes();
   }

   void RadialTruncation::correctGridSizes()
   {
      if(SimulationConstants::isCSCSGrid())
      {
         this->mNr += 2;
      }
   }
}
