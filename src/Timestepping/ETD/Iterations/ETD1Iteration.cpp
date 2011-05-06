/** \file ETD1Iteration.cpp
 *  \brief Implemenation of the simple ETD1 iteration step. This is also the A step for ETD2RK method. 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/Iterations/ETD1Iteration.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD1Iteration::ETD1Iteration(EPMFloat nFactor, SmartETDOperators pOpM0, SmartETDOperators pOpM1)
      : mNFactor(nFactor), mpOpM0(pOpM0), mpOpM1(pOpM1)
   {
   }

   void ETD1Iteration::compute(ETD1Iteration::ScalarType &rVar, ETD1Iteration::ScalarType &rNTerms)
   {
      // Get number of harmonic degrees
      int nL = rVar.trunc()->local()->spec()->nL();

      // Get minimal degree index (not l=0)
      int l0 = rVar.minL();

      // loop over degrees 
      for(int l = l0; l < nL; ++l)
      {  
         rNTerms.rLShell(l) *= this->mNFactor;
         this->mpOpM1->multiplyOrders(rNTerms.rLShell(l), rNTerms.lshell(l),l);
         this->mpOpM0->affineOrders(rVar.rLShell(l), rVar.lshell(l), rNTerms.lshell(l), l);
         this->mpOpM0->extendOrders(rVar.rLShell(l), l);
      }
   }

}
