/** \file ThetaPredictor.cpp
 *  \brief Implemenation of the predictor step of the theta method
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/PredictorCorrector/Theta/ThetaPredictor.hpp"

// Project includes
//

namespace EPMDynamo {

   ThetaPredictor::ThetaPredictor(SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS)
      : mpOpLHS(pOpLHS), mpOpRHS(pOpRHS)
   {
   }

   void ThetaPredictor::compute(ThetaPredictor::ScalarType &rVar, ThetaPredictor::ScalarType &rNTerms)
   {
      // Set the RHS part of predictor step
      this->setRHS(rVar, rVar, rNTerms);

      // Solve predictor equations
      this->solve(rVar);
   }

   void ThetaPredictor::setRHS(ThetaPredictor::ScalarType& rVar, const ThetaPredictor::ScalarType& oldVar, const ThetaPredictor::ScalarType& nTerms)
   {
      // Get number of harmonic degrees
      int nL = this->mpOpRHS->trunc()->local()->spec()->nL();
      const int l0 = rVar.minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Compute RHS part of predictor timestep equation for given degree
         this->mpOpRHS->affineOrders(rVar.rLShell(l), oldVar.lshell(l), nTerms.lshell(l), l);
      }
   }

   void ThetaPredictor::solve(ThetaPredictor::ScalarType& rVar)
   {
      // Solve for unknown variable
      this->mpOpLHS->solve(rVar);
   }
}
