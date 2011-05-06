/** \file ETD2RKMethod.cpp
 *  \brief Implementation of the ETD2RK method (without influence matrix)
 */


// Configuration includes
//

// System includes
//

// External includes
//

// Class includes
//
#include "Timestepping/ETD/ETD2RK/ETD2RKMethod.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD2RKMethod::ETD2RKMethod(EPMFloat a, EPMFloat b, const ETD2RKMethod::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : IterativeSchemeBase(a, b, basis, tsteps, pTrunc, hasL0), mETD2(b/a, pTrunc, hasL0)
   {
   }

   void ETD2RKMethod::addBC(SmartBC pBC)
   {
      // add boundary condition to operators
      this->mETD2.addBC(pBC);
   }

   void ETD2RKMethod::init()
   {
      // initialise pointers
      this->initStorage();

      // Initialise the operators
      this->initOperators();
   }

   void ETD2RKMethod::initOperators()
   {
      // initialise the operators
      this->mETD2.initOperators(this->mrBasis);
   }

   void ETD2RKMethod::updateTimeMatrices()
   {
      // Update the time depended matrices
      this->mETD2.update(this->rTSParams().dt(), this->mrBasis);
   }

   void ETD2RKMethod::initStorage()
   {
      // Create intermediate value a computation step
      EPMSHARED_PTR<ETD1Iteration > pItA(new ETD1Iteration (1.0/this->mA, this->mETD2.pEtdF(0), this->mETD2.pEtdF(1)));

      // Create timestep computation step
      EPMSHARED_PTR<ETD2RKTimestep > pItTimestep(new ETD2RKTimestep (1.0/this->mA, this->pOldNTerms(), this->mETD2.pEtdF(2)));

      // Add required ETD steps
         // Add intermediate value A computation
      this->mSteps.push_back(pItA);
         // Add timestep step
      this->mSteps.push_back(pItTimestep);
   }

}
