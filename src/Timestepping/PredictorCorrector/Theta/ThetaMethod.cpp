/** \file ThetaMethod.cpp
 *  \brief Implementation of the predictor-corrector theta method (without influence matrix)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/PredictorCorrector/Theta/ThetaMethod.hpp"

// Project includes
//

namespace EPMDynamo {

   ThetaMethod::ThetaMethod(EPMFloat a, EPMFloat b, const ThetaMethod::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : IterativeSchemeBase(a, b, basis, tsteps, pTrunc, hasL0), mpLHS(new ThetaLHSTOperatorSet(a, b, basis, pTrunc, hasL0)), mpRHS(new ThetaRHSTOperatorSet(a, b, basis, pTrunc, hasL0))
   {
   }

   void ThetaMethod::addBC(SmartBC pBC)
   {
      // Add boundary condition to LHS of theta method
      this->mpLHS->addBC(pBC);

      // Add boundary condition to RHS of theta method
      this->mpRHS->addBC(pBC);
   }

   void ThetaMethod::init()
   {
      // initialise pointers
      this->initStorage();

      // initialise the operators
      this->initOperators();
   }

   void ThetaMethod::initOperators()
   {
      // Create LHS operators and LU factorise them
      this->mpLHS->initOperators();

      // Create RHS operators
      this->mpRHS->initOperators();
   }

   void ThetaMethod::updateTimeMatrices()
   {
      // Updated the LHS operator due to new timestep
      this->mpLHS->update(this->rTSParams().dt());

      // Updated the RHS operator due to new timestep
      this->mpRHS->update(this->rTSParams().dt());
   }

   void ThetaMethod::initStorage()
   {
      //
      // Create intermediate storage
      EPMSHARED_PTR<ScalarType> pPrevious(new ScalarType (this->mpLHS->trunc(), this->mpLHS->hasL0()));
      this->mTmpScalars.push_back(pPrevious);

      // Create the predictor computation step
      EPMSHARED_PTR<ThetaPredictor> pItP(new ThetaPredictor (this->mpLHS, this->mpRHS));

      // Create the corrector computation step
      EPMSHARED_PTR<ThetaCorrector> pItC(new ThetaCorrector (pPrevious, this->mpLHS, this->mpRHS));

      //
      // Add the required computation steps
      
      // Add the predictor computation
      this->mSteps.push_back(pItP);
      // Add a corrector step
      this->mSteps.push_back(pItC);
   }

   void ThetaMethod::doIteration(ThetaMethod::ScalarType& rVar, ThetaMethod::ScalarType& rNTerms)
   {
      // Do normal iteration
      IterativeSchemeBase::doIteration(rVar, rNTerms);

      // Update values stored from previous iteration
      this->updatePrevious(rNTerms);
   }

}
