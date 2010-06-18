/** \file PCTimestepConfig.cpp
 *  \brief Source of the specialised timestep config for a predictor/corrector scheme
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/PredictorCorrector/PCTimestepConfig.hpp"
#include "Timestepping/TimestepConfig.hpp"

// Project includes
//

namespace EPMDynamo {

   const DynamoFloat PCTimestepConfig::MAX_CORRECTOR_TSTEP = 0.5;

   const DynamoFloat PCTimestepConfig::MAX_CORRECTION_NORM = TimestepConfig::TIMESTEP_ERROR_EPSILON;

   PCTimestepConfig::PCTimestepConfig()
   {
   }

}
