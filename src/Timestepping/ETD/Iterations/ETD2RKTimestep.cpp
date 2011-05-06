/** \file ETD2RKTimestep.cpp
 *  \brief Implementation of the ETD2RK timestep step
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/Iterations/ETD2RKTimestep.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD2RKTimestep::ETD2RKTimestep(EPMFloat nFactor, ETD2RKTimestep::SmartScalarType pOldN, ETD2RKTimestep::SmartETDOperators pOpM2)
      : mNFactor(nFactor), mpOldN(pOldN), mpOpM2(pOpM2)
   {
   }

   void ETD2RKTimestep::prepare(ETD2RKTimestep::ScalarType &rVar, ETD2RKTimestep::ScalarType &rNTerms)
   {
      // Get number of harmonic degrees
      int nL = rVar.trunc()->local()->spec()->nL();

      // Get minimal degree index (not l=0)
      int l0 = rVar.minL();

      // loop over degrees 
      for(int l = l0; l < nL; ++l)
      {
         rNTerms.rLShell(l) -= this->mpOldN->lshell(l);
         rNTerms.rLShell(l) *= this->mNFactor;
      }
   }

   void ETD2RKTimestep::compute(ETD2RKTimestep::ScalarType &rVar, ETD2RKTimestep::ScalarType &rNTerms)
   {
      // Get number of harmonic degrees
      int nL = rVar.trunc()->local()->spec()->nL();

      // Get minimal degree index (not l=0)
      int l0 = rVar.minL();

      // loop over degrees 
      for(int l = l0; l < nL; ++l)
      {
         this->mpOpM2->multiplyOrders(rNTerms.rLShell(l), rNTerms.lshell(l), l);
         rVar.rLShell(l) += rNTerms.lshell(l);
         this->mpOpM2->extendOrders(rVar.rLShell(l), l);
      }
   }

}
