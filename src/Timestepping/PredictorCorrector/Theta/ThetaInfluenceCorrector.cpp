/** \file ThetaInfluenceCorrector.hpp
 *  \brief Implemenation of the corrector step of the theta method
 */

#ifndef THETAINFLUENCECORRECTOR_HPP
#define THETAINFLUENCECORRECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/PredictorCorrector/Theta/ThetaInfluenceCorrector.hpp"

// Project includes
//

namespace EPMDynamo {

   ThetaInfluenceCorrector::ThetaInfluenceCorrector(ThetaInfluenceCorrector::SmartInfluenceMatrix pInfluence, ThetaInfluenceCorrector::SmartScalarType pPrevious, ThetaInfluenceCorrector::SmartLHSOperators pOpLHS, ThetaInfluenceCorrector::SmartRHSOperators pOpRHS)
      : ThetaCorrector(pPrevious, pOpLHS, pOpRHS), mpInfluence(pInfluence)
   {
   }

   void ThetaInfluenceCorrector::compute(ThetaInfluenceCorrector::ScalarType &rVar, ThetaInfluenceCorrector::ScalarType &rNTerms)
   {
      // Set the RHS part of corrector step
      this->setRHS(rNTerms);

      // Solve corrector equations
      this->solve(rNTerms);

      // Use corrector solution
      this->useCorrection(rVar, rNTerms);
   }

   void ThetaInfluenceCorrector::setRHS(ThetaInfluenceCorrector::ScalarType& rNTerms)
   {
      // Solve influenc matrix equation
      this->mpInfluence->solve(rNTerms);
      
      // Compute RHS part of corrector timestep equation
      ThetaCorrector::setRHS(rNTerms);
   }

   void ThetaInfluenceCorrector::solve(ThetaInfluenceCorrector::ScalarType& rNTerms)
   {
      // Get the correction to the unknown variable
      ThetaCorrector::solve(rNTerms);

      // Include kernel influence
      this->mpInfluence->addKernel(rNTerms);
   }
}
