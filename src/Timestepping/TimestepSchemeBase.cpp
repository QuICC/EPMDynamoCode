/** \file TimestepSchemeBase.cpp
 *  \brief Base for any implementation of timestep scheme
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/TimestepSchemeBase.hpp"

// Project includes
//

namespace EPMDynamo {

   TimestepSchemeBase::TimestepSchemeBase(TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : mrTStepParams(tsteps), mpOldVar(new ScalarType(pTrunc, hasL0)), mpOldNTerms(new ScalarType(pTrunc, hasL0)) 
   {
   }

   void TimestepSchemeBase::storeOld(const TimestepSchemeBase::ScalarType& var, const TimestepSchemeBase::ScalarType& nTerms)
   {
      // Get number of harmonic degrees
      int nL = this->oldNTerms().nL();
      // Get minimum harmonic degrees
      const int l0 = this->oldNTerms().minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Store previous timestep variables
         this->rOldVar().rLShell(l) = var.lshell(l);

         // Store previous timestep nterms
         this->rOldNTerms().rLShell(l) = nTerms.lshell(l);
      }
   }

   void TimestepSchemeBase::restoreOld(TimestepSchemeBase::ScalarType& rVar, TimestepSchemeBase::ScalarType& nTerms)
   {
      // Get number of harmonic degrees
      int nL = this->oldNTerms().nL();
      // Get minimu harmonic degrees
      const int l0 = this->oldNTerms().minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Recover previous timestep variables
         rVar.rLShell(l) = this->oldVar().lshell(l);

         // Recover previous timestep nterms
         nTerms.rLShell(l) = this->oldNTerms().lshell(l);
      }
   }
}
