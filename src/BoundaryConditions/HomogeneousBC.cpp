/** \file HomogeneousBC.cpp
 *  \brief Source of the Homogeneous boundary condition implemenation
 */

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/HomogeneousBC.hpp"

// Project includes
//

namespace EPMDynamo {

   HomogeneousBC::HomogeneousBC(const int nL, const int nX)
      : BoundaryCondition(nL, nX)
   {
      // These are all homogeneous boundary conditions
      this->mIsHomogeneous = true;
   }
}
