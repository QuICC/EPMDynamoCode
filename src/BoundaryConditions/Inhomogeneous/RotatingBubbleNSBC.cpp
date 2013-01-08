/** \file RotatingBubbleNSBC.cpp
 *  \brief Implementation of the basis first rotating bubble poloidal boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/Inhomogeneous/RotatingBubbleNSBC.hpp"

// Project includes
//

namespace EPMDynamo {

   RotatingBubbleNSBC::RotatingBubbleNSBC(const RotatingBubbleNSBC::BasisType &basis)
      : BoundaryCondition(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      fillLHSBCValues(basis);
   }

   void RotatingBubbleNSBC::fillLHSBCValues(const RotatingBubbleNSBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

   EPMFloat RotatingBubbleNSBC::getRealRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return 0.0;
      } else
      {
         return 0.0;
      }
   }

   EPMFloat RotatingBubbleNSBC::getImagRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return 0.0;
      } else
      {
         return 0.0;
      }
   }

   EPMComplex RotatingBubbleNSBC::getRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return EPMComplex(0.0, 0.0);
      } else
      {
         return EPMComplex(0.0);
      }
   }

}
