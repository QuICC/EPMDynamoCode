/** \file DRadialBC.cpp
 *  \brief Implementation of the homogeneous first radial derivative boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"

// Project includes
//

namespace EPMDynamo {

   DRadialBC::DRadialBC(const DRadialBC::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void DRadialBC::fillLHSBCValues(const DRadialBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bdiff(1);
      }
   }
}
