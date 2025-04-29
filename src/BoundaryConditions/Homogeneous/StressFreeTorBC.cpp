/** \file StressFreeTorBC.cpp
 *  \brief Implementation of the stress free velocity boundary condition on the toroidal scalar
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

// Project includes
//

namespace EPMDynamo {

   StressFreeTorBC::StressFreeTorBC(const StressFreeTorBC::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void StressFreeTorBC::fillLHSBCValues(const StressFreeTorBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bdiff(1) - basis.at(l).bpoly();
      }
   }

}
