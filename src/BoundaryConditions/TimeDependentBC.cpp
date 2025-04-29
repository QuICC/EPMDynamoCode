/** \file TimeDependentBC.cpp
 *  \brief Source of the TimeDependent boundary condition implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/TimeDependentBC.hpp"

// Project includes
//

namespace EPMDynamo {

   TimeDependentBC::TimeDependentBC(const int nL, const int nX, const TimestepParameters &tsParams)
      : BoundaryCondition(nL, nX), mrTSParams(tsParams)
   {
   }
}
