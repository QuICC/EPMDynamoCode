/** \file ThetaCorrector.cpp
 *  \brief Implemenation of the corrector step of the theta method
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/PredictorCorrector/Theta/ThetaCorrector.hpp"

// Project includes
//

namespace EPMDynamo {

   ThetaCorrector::ThetaCorrector(ThetaCorrector::SmartScalarType pPrevious, ThetaCorrector::SmartLHSOperators pOpLHS, ThetaCorrector::SmartRHSOperators pOpRHS)
      : mpPreviousNTerms(pPrevious), mpOpLHS(pOpLHS), mpOpRHS(pOpRHS)
   {
   }

   void ThetaCorrector::compute(ThetaCorrector::ScalarType &rVar, ThetaCorrector::ScalarType &rNTerms)
   {
      // Set the RHS part of corrector step
      this->setRHS(rNTerms);

      // Solve corrector equations
      this->solve(rNTerms);

      // Use corrector solution
      this->useCorrection(rVar, rNTerms);
   }

   void ThetaCorrector::setRHS(ThetaCorrector::ScalarType& rNTerms)
   {
      // Get number of harmonic degrees
      int nL = this->mpOpRHS->trunc()->local()->spec()->nL();
      const int l0 = this->mpPreviousNTerms->minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         rNTerms.rLShell(l) = ThetaTraits::theta*(rNTerms.lshell(l) - this->mpPreviousNTerms->lshell(l));
      }
   }

   void ThetaCorrector::solve(ThetaCorrector::ScalarType& rNTerms)
   {
      // Get the correction to the unknown variable
      this->mpOpLHS->solveZero(rNTerms);
   }

   void ThetaCorrector::useCorrection(ThetaCorrector::ScalarType& rVar, const ThetaCorrector::ScalarType& nTerms)
   {
      // Add correction to unknown
      this->addCorrection(rVar, nTerms);
   }

   void ThetaCorrector::addCorrection(ThetaCorrector::ScalarType& rVar, const ThetaCorrector::ScalarType& corr) const
   {
      // Get number of harmonic degrees and minimal index
      int nL = rVar.nL();
      int l0 = rVar.minL();

      // loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Add corrector correction to solution
         rVar.rLShell(l) += corr.lshell(l);
      }
   }
}
