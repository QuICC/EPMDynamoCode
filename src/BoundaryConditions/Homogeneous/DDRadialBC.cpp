/** \file DDRadialBC.cpp
 *  \brief Implementation of the homogeneous second radial derivative boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"

// Project includes
//

namespace EPMDynamo {

   DDRadialBC::DDRadialBC(const DDRadialBC::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void DDRadialBC::fillLHSBCValues(const DDRadialBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bdiff(2);
      }
   }
}
