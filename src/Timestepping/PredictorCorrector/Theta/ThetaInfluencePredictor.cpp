/** \file ThetaInfluencePredictor.cpp
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
#include "Timestepping/PredictorCorrector/Theta/ThetaInfluencePredictor.hpp"

// Project includes
//

namespace EPMDynamo {

   ThetaInfluencePredictor::ThetaInfluencePredictor(SmartInfluenceMatrix pInfluence, SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS)
      : ThetaPredictor(pOpLHS, pOpRHS), mpInfluence(pInfluence)
   {
   }

   void ThetaInfluencePredictor::compute(ThetaInfluencePredictor::ScalarType &rVar, ThetaInfluencePredictor::ScalarType &rNTerms)
   {
      // Set the RHS part of predictor step
      this->setRHS(rVar, rVar, rNTerms);

      // Solve predictor equations
      this->solve(rVar);
   }

   void ThetaInfluencePredictor::setRHS(ThetaInfluencePredictor::ScalarType& rVar, const ThetaInfluencePredictor::ScalarType& oldVar, ThetaInfluencePredictor::ScalarType& rNTerms)
   {
      // Solve influence matrix equation
      this->mpInfluence->solve(rNTerms);

      // Compute RHS part of predictor timestep equation
      ThetaPredictor::setRHS(rVar, oldVar, rNTerms);
   }

   void ThetaInfluencePredictor::solve(ThetaInfluencePredictor::ScalarType& rVar)
   {
      // Solve for unknown variable
      ThetaPredictor::solve(rVar);

      // Include kernel influence
      this->mpInfluence->addKernel(rVar);
   }
}
