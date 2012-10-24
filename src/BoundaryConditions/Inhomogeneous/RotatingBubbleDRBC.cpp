/** \file RotatingBubbleDRBC.cpp
 *  \brief Implementation of the basis second rotating bubble poloidal boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/Inhomogeneous/RotatingBubbleDRBC.hpp"

// Project includes
//

namespace EPMDynamo {

   RotatingBubbleDRBC::RotatingBubbleDRBC(const RotatingBubbleDRBC::BasisType &basis)
      : BoundaryCondition(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void RotatingBubbleDRBC::fillLHSBCValues(const RotatingBubbleDRBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bdiff(1);
      }
   }

   EPMFloat RotatingBubbleDRBC::getRealRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return 0.69099*1.0364;
      } else
      {
         return 0.0;
      }
   }

   EPMFloat RotatingBubbleDRBC::getImagRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return 0.0;
      } else
      {
         return 0.0;
      }
   }

   EPMComplex RotatingBubbleDRBC::getRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return EPMComplex(0.69099*1.0364, 0.0);
      } else
      {
         return EPMComplex(0.0);
      }
   }
}
