/** \file TimestepConfig.cpp
 *  \brief Source of some timestepping constants
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/TimestepConfig.hpp"

// Project includes
//

namespace EPMDynamo {

   const EPMFloat TimestepConfig::MIN_TIMESTEP = 1e-10;

   EPMFloat TimestepConfig::MAXIMUM_TIMESTEP = 1.0;

   const EPMFloat TimestepConfig::COURANT_NUMBER = 0.8;

   const EPMFloat TimestepConfig::TIMESTEP_MARGIN = 0.03;

   const EPMFloat TimestepConfig::TIMESTEP_ERROR_THETA = 0.5;

   const EPMFloat TimestepConfig::TIMESTEP_ERROR_TOLERANCE = 1e-8;

   const EPMFloat TimestepConfig::TIMESTEP_ERROR_EPSILON = TimestepConfig::TIMESTEP_ERROR_THETA*TimestepConfig::TIMESTEP_ERROR_TOLERANCE;

   const EPMFloat TimestepConfig::TIMESTEP_RELERROR_SCALING = 1.0;
   // const EPMFloat TimestepConfig::TIMESTEP_RELERROR_SCALING = 0.0;
   
   const EPMFloat TimestepConfig::TIMESTEP_MAX_RATIO = 1.602;

   void TimestepConfig::overruleMaxTimestep(EPMFloat newMax)
   {
      // Set new maximum timestep
      MAXIMUM_TIMESTEP = newMax;

      // Print message about change of maximum timestep
      std::cout << "MAXIMUM TIMESTEP has been overruled and set to :" << MAXIMUM_TIMESTEP << std::endl;
      std::cout << std::endl;
   }

}
