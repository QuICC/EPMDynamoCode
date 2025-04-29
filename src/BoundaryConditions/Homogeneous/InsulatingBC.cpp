/** \file InsulatingBC.cpp
 *  \brief Implementation of the magnetic insulating boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/Homogeneous/InsulatingBC.hpp"

// Project includes
//

namespace EPMDynamo {

   InsulatingBC::InsulatingBC(const InsulatingBC::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void InsulatingBC::fillLHSBCValues(const InsulatingBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bdiff(1) + static_cast<EPMFloat>(basis.at(l).l()+1)*basis.at(l).bpoly();
      }
   }

}
