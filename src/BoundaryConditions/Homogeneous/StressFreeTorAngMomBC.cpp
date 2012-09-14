/** \file StressFreeTorAngMomBC.cpp
 *  \brief Implementation of the stress free velocity boundary condition on the toroidal scalar include angular momentum conservation for l=1
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "BoundaryConditions/Homogeneous/StressFreeTorAngMomBC.hpp"

// Project includes
//

namespace EPMDynamo {

   StressFreeTorAngMomBC::StressFreeTorAngMomBC(const StressFreeTorAngMomBC::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void StressFreeTorAngMomBC::fillLHSBCValues(const StressFreeTorAngMomBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = basis.at(l).bdiff(1) - basis.at(l).bpoly();
         if(basis.at(l).l() == 1)
         {
            this->rLHSBCValues(l).array() *= 0.25;

            for(int n = 0; n < basis.polyN(); ++n)
            {
               this->rLHSBCValues(l)(n) += basis.at(l).eWeights()(0,n);
            }
         }
         std::cerr << this->rLHSBCValues(l).transpose() << std::endl;
      }
   }

}
