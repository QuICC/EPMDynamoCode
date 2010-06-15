/** \file SimulationConstants.cpp
 *  \brief Source of the simulation wide static constants class
 */

// System includes
//

// External includes
//

// Class include
//
#include "Simulations/SimulationConstants.hpp"

// Project includes
//
#include "General/EPMException.hpp"

namespace EPMDynamo {

   bool SimulationConstants::msIsCSCS = false;

   bool SimulationConstants::msIsGridLocked = false;

   void SimulationConstants::setCSCSGrid()
   {
      if(SimulationConstants::isGridLocked())
      {
         throw EPMException("SimulationConstants::setCSCSGrid", "Tried to change grid type while it has already been used");
      } else
      {
         SimulationConstants::msIsCSCS = true;
      }
   }

   void SimulationConstants::lockGridType()
   {
      SimulationConstants::msIsGridLocked = true;
   }

}
