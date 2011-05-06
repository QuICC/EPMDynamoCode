/** \file ThetaInfluenceMethod.cpp
 *  \brief Implementation of the predictor-corrector theta method with influence matrix
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/PredictorCorrector/Theta/ThetaInfluenceMethod.hpp"

// Project includes
//
#include "General/EPMException.hpp"

namespace EPMDynamo {

   ThetaInfluenceMethod::ThetaInfluenceMethod(EPMFloat a, EPMFloat b, const ThetaInfluenceMethod::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : ThetaMethod(a, b, basis, tsteps, pTrunc, hasL0), mInfluenceNBC(-2), mpInfluence(new InfluenceMatrix (pTrunc, basis, hasL0))
   {
   }

   void ThetaInfluenceMethod::addBC(SmartBC pBC)
   {
      // Propagate only the first boundary condition to standard implementation
      if(this->mInfluenceNBC == -2)
      {
         ThetaMethod::addBC(pBC);
      }

      // Add boundary condition to influence matrix
      this->mpInfluence->addBC(pBC);

      // Increment number of implemented influence matrix BCs
      ++this->mInfluenceNBC;
   }

   void ThetaInfluenceMethod::init()
   {
      // initialise pointers
      this->initStorage();

      // Initialise the operators
      this->initOperators();

      // Initialise the influence matrix
      this->initInfluence();
   }

   void ThetaInfluenceMethod::updateTimeMatrices()
   {
      // Update the timestep operators
      ThetaMethod::updateTimeMatrices();

      // Update the influence matrix solution
      this->updateInfluence();
   }

   void ThetaInfluenceMethod::initStorage()
   {
      //
      // Create intermediate storage
      EPMSHARED_PTR<ScalarType > pPrevious(new ScalarType (this->mpLHS->trunc(), this->mpLHS->hasL0()));
      this->mTmpScalars.push_back(pPrevious);

      // Create the predictor computation step
      EPMSHARED_PTR<ThetaInfluencePredictor > pItP(new ThetaInfluencePredictor (this->mpInfluence, this->mpLHS, this->mpRHS));

      // Create the corrector computation step
      EPMSHARED_PTR<ThetaInfluenceCorrector > pItC(new ThetaInfluenceCorrector (this->mpInfluence, pPrevious, this->mpLHS, this->mpRHS));

      //
      // Add the required computation steps
      
      // Add the predictor computation
      this->mSteps.push_back(pItP);
      // Add a corrector step
      this->mSteps.push_back(pItC);
   }

   void ThetaInfluenceMethod::initInfluence()
   {
      if(this->mInfluenceNBC  == 0)
      {
         // Initialise the operators
         this->mpInfluence->initOperators();

         // Compute the operators
         this->mpInfluence->computeOperators();
      } else
      {
         throw EPMException("ThetaInfluenceMethod::initInfluence", "Tried to initialise with wrong number of BCs");
      }
   }

   void ThetaInfluenceMethod::updateInfluence()
   {
      // Get size of radial truncation
      int nN = this->mpLHS->trunc()->sim()->rad()->nN();
      // Get number of harmonic degrees
      int nL = this->mpLHS->trunc()->local()->spec()->nL();

      // Get minimal degree index (not l=0)
      int l0 = ! this->mpLHS->trunc()->local()->spec()->lArray()(0);

      // Create temporary storage
      Array tmp(nN);

      // loop over degrees 
      for(int l = l0; l < nL; ++l)
      {
         // Initialise influence matrix solution to kernel r^l
         tmp.setConstant(0.0);
         tmp(0) = 1.0;

         // Compute the kernel influence
         this->mpLHS->solveZeroVector(tmp, l);

         // Store solution from influence matrix
         this->mpInfluence->storeKernelBC(tmp, l);
      }
   }

}
