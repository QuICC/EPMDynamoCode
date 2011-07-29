/** \file LibrationBC.cpp
 *  \brief Implementation of the time dependent poloidal libration boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/TimeDependent/LibrationBC.hpp"

// Project includes
//

namespace EPMDynamo {

   LibrationBC::LibrationBC(const LibrationBC::BasisType &basis, const TimestepParameters &tsParams)
      : TimeDependentBC(basis.basisN(), basis.polyN(), tsParams)
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   EPMFloat LibrationBC::evolvingRealFactor() const
   {
      EPMFloat epsilon = 0.001;
      EPMFloat f = 0.99;

      return epsilon*std::cos(f*this->mrTSParams.time());
   }

   EPMFloat LibrationBC::evolvingImagFactor() const
   {
      return 0.0;
   }

   void LibrationBC::fillLHSBCValues(const LibrationBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

   EPMFloat LibrationBC::getRealRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return -0.5*this->evolvingRealFactor();
      } else
      {
         return 0.0;
      }
   }

   EPMFloat LibrationBC::getImagRHSBC(const int l, const int m) const
   {
      return 0.0;
   }

   EPMComplex LibrationBC::getRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return EPMComplex(-0.5*this->evolvingRealFactor(), 0.0);
      } else
      {
         return EPMComplex(0.0);
      }
   }
}
