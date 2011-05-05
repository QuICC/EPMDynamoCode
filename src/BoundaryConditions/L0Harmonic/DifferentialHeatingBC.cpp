/** \file DifferentialHeatingBC.cpp
 *  \brief Implementation of the l=0 differential heating boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "BoundaryConditions/L0Harmonic/DifferentialHeatingBC.hpp"

// Project includes
//

namespace EPMDynamo {

   DifferentialHeatingBC::DifferentialHeatingBC(const EPMFloat bcValue, const DifferentialHeatingBC::BasisType &basis)
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

   void DifferentialHeatingBC::fillLHSBCValues(const DifferentialHeatingBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

   void DifferentialHeatingBC::setBCValue(const EPMFloat bcValue)
   {
      this->mValue = bcValue;
   }

}
