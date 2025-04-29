/** \file PrecessionFrameBC.cpp
 *  \brief Implementation of changed frame of referecen toroidal precession boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/Inhomogeneous/PrecessionFrameBC.hpp"

// Project includes
//

namespace EPMDynamo {

   PrecessionFrameBC::PrecessionFrameBC(const PrecessionFrameBC::BasisType &basis)
      : BoundaryCondition(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void PrecessionFrameBC::fillLHSBCValues(const PrecessionFrameBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

   EPMFloat PrecessionFrameBC::getRealRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 0)
      {
         return 1.0;
      } else
      {
         return 0.0;
      }
   }

   EPMFloat PrecessionFrameBC::getImagRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 0)
      {
         return 0.0;
      } else
      {
         return 0.0;
      }
   }

   EPMComplex PrecessionFrameBC::getRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 0)
      {
         return EPMComplex(1.0, 0.0);
      } else
      {
         return EPMComplex(0.0);
      }
   }
}
