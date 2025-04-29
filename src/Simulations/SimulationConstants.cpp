/** \file SimulationConstants.cpp
 *  \brief Source of the simulation wide static constants class
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "Simulations/SimulationConstants.hpp"

// Project includes
//

namespace EPMDynamo {

   bool SimulationConstants::msIsCSCS = false;

   bool SimulationConstants::msIsGridLocked = false;

   void SimulationConstants::setCSCSGrid()
   {
      assert(! SimulationConstants::isGridLocked());

      SimulationConstants::msIsCSCS = true;
   }

   void SimulationConstants::lockGridType()
   {
      SimulationConstants::msIsGridLocked = true;
   }

}
