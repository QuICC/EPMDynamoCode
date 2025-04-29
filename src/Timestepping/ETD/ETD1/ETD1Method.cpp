/** \file ETD1Method.cpp
 *  \brief Implementation of the ETD1 method (without influence matrix)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/ETD1/ETD1Method.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD1Method::ETD1Method(EPMFloat a, EPMFloat b, const ETD1Method::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : IterativeSchemeBase(a, b, basis, tsteps, pTrunc, hasL0), mETD1(b/a, pTrunc, hasL0)
   {
   }

   void ETD1Method::addBC(SmartBC pBC)
   {
      // add boundary condition to operators
      this->mETD1.addBC(pBC);
   }

   void ETD1Method::init()
   {
      // initialise pointers
      this->initStorage();

      // Initialise the operators
      this->initOperators();
   }

   void ETD1Method::initOperators()
   {
      // initialise the operators
      this->mETD1.initOperators(this->mrBasis);
   }

   void ETD1Method::updateTimeMatrices()
   {
      // Update the time depended matrices
      this->mETD1.update(this->rTSParams().dt(), this->mrBasis);
   }

   void ETD1Method::initStorage()
   {
      // Create intermediate value a computation step
      EPMSHARED_PTR<ETD1Iteration > pIt(new ETD1Iteration (1.0/this->mA, this->mETD1.pEtdF(0), this->mETD1.pEtdF(1)));

      // Add required ETD steps
         // Add intermediate value A computation
      this->mSteps.push_back(pIt);
   }

}
