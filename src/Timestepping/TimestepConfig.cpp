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

   const DynamoFloat TimestepConfig::MIN_TIMESTEP = 1e-10;

   DynamoFloat TimestepConfig::MAXIMUM_TIMESTEP = 1.0;

   const DynamoFloat TimestepConfig::COURANT_NUMBER = 0.8;

   const DynamoFloat TimestepConfig::TIMESTEP_MARGIN = 0.03;

   const DynamoFloat TimestepConfig::TIMESTEP_ERROR_THETA = 0.5;

   const DynamoFloat TimestepConfig::TIMESTEP_ERROR_TOLERANCE = 1e-2;

   const DynamoFloat TimestepConfig::TIMESTEP_ERROR_EPSILON = TimestepConfig::TIMESTEP_ERROR_THETA*TimestepConfig::TIMESTEP_ERROR_TOLERANCE;

   TimestepConfig::TimestepConfig()
   {
   }

   void TimestepConfig::overruleMaxTimestep(DynamoFloat newMax)
   {
      MAXIMUM_TIMESTEP = newMax;

      std::cout << "MAXIMUM TIMESTEP has been overruled and set to :" << MAXIMUM_TIMESTEP << std::endl;
   }

}
