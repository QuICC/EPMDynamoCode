/** \file IterativeSchemeBase.cpp
 *  \brief Base of the implementation of an iterative timestep scheme
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class includes
//
#include "Timestepping/IterativeSchemeBase.hpp"

// Project includes
//
#include "Timestepping/ErrorControl.hpp"

namespace EPMDynamo {

   IterativeSchemeBase::IterativeSchemeBase(EPMFloat a, EPMFloat b, const IterativeSchemeBase::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : TimestepSchemeBase(tsteps, pTrunc, hasL0), mA(a), mB(b), mrBasis(basis)
   {
   }

   void IterativeSchemeBase::resetIterations()
   {
      this->mCurrentIt = this->mSteps.begin();
   }

   void IterativeSchemeBase::doIteration(IterativeSchemeBase::ScalarType& rVar, IterativeSchemeBase::ScalarType& rNTerms)
   {
      // Do step computation
      (*this->mCurrentIt)->compute(rVar, rNTerms);

      // Update the error from timestep
      if((*this->mCurrentIt)->providesError())
      {
         this->rTSParams().updateError(ErrorControl::errorNorm(rNTerms, this->oldVar(), this->rTSParams().error()));
      }

      // Go one iteration further (this allows to stay on the same iteration several times)
      if((*this->mCurrentIt)->doNextIteration())
      {
         // Go forward one step
         ++this->mCurrentIt;
      }
   }
}
