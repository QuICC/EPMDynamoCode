/** \file L0HarmonicBC.cpp
 *  \brief Source of the l=0 harmonic boundary condition implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/L0HarmonicBC.hpp"

// Project includes
//

namespace EPMDynamo {

   L0HarmonicBC::L0HarmonicBC(const int nL, const int nX, const bool hasL0)
      : BoundaryCondition(nL, nX), mHasL0(hasL0), mValue(0.0)
   {
   }
}
