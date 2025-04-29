/** \file ZeroBC.cpp
 *  \brief Implementation of the basis Zero value boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"

// Project includes
//

namespace EPMDynamo {

   ZeroBC::ZeroBC(const ZeroBC::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      fillLHSBCValues(basis);
   }

   void ZeroBC::fillLHSBCValues(const ZeroBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

}
