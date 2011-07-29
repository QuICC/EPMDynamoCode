/** \file PrecessionBC.cpp
 *  \brief Implementation of the time dependent toroidal precession boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/TimeDependent/PrecessionBC.hpp"

// Project includes
//

namespace EPMDynamo {

   PrecessionBC::PrecessionBC(const PrecessionBC::BasisType &basis, const TimestepParameters &tsParams)
      : TimeDependentBC(basis.basisN(), basis.polyN(), tsParams)
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   EPMFloat PrecessionBC::evolvingRealFactor() const
   {
      return 0.0;
   }

   EPMFloat PrecessionBC::evolvingImagFactor() const
   {
      return 0.0;
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
      return 0.0;
   }

   EPMComplex PrecessionBC::getRHSBC(const int l, const int m) const
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
