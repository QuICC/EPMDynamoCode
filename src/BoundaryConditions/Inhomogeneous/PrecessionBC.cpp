/** \file PrecessionBC.cpp
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
#include "BoundaryConditions/Inhomogeneous/PrecessionBC.hpp"

// Project includes
//

namespace EPMDynamo {

   PrecessionBC::PrecessionBC(const PrecessionBC::BasisType &basis)
      : BoundaryCondition(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void PrecessionBC::fillLHSBCValues(const PrecessionBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

   EPMFloat PrecessionBC::getRealRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 0)
      {
         return 1.0;
      } else
      {
         return 0.0;
      }
   }

   EPMFloat PrecessionBC::getImagRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 0)
      {
         return 1.0;
      } else
      {
         return 0.0;
      }
   }

   EPMComplex PrecessionBC::getRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 0)
      {
         return EPMComplex(1.0);
      } else
      {
         return EPMComplex(0.0);
      }
   }
}
