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

   const EPMFloat PCTimestepConfig::MAX_CORRECTION_NORM = TimestepConfig::TIMESTEP_ERROR_EPSILON;

   const EPMFloat PCTimestepConfig::MAX_CORRECTOR_ERROR = 5*MAX_CORRECTION_NORM;

}
