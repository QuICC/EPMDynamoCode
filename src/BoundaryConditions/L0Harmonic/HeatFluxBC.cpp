/** \file HeatFluxBC.cpp
 *  \brief Implementation of the l=0 constant heat flux boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "BoundaryConditions/L0Harmonic/HeatFluxBC.hpp"

// Project includes
//

namespace EPMDynamo {

   HeatFluxBC::HeatFluxBC(const EPMFloat bcValue, const HeatFluxBC::BasisType &basis)
      : L0HarmonicBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);

      // Set boundary value
      this->setBCValue(bcValue);

      // Check if boundary condition is homogeneous
      if(bcValue == 0.0)
      {
         this->mIsHomogeneous = true;
      }
   }

   void HeatFluxBC::fillLHSBCValues(const HeatFluxBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bdiff(1);
      }
   }

   void HeatFluxBC::setBCValue(const EPMFloat bcValue)
   {
      this->mValue = bcValue;
   }

}
