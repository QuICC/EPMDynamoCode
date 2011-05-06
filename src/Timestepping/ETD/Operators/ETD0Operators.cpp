/** \file ETD0Operators.cpp
 *  \brief Implementation of the ETD0 operators
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/Operators/ETD0Operators.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD0Operators::ETD0Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : ETDNOperators<1>(c, pTrunc, hasL0)
   {
   }

   void ETD0Operators::update(const EPMFloat dt, const ETD0Operators::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   void ETD0Operators::createOperators(const EPMFloat h, const ETD0Operators::BasisType &basis)
   {
      // Update the required scaling power
      this->updateScalings(h);

      // Loop over all degrees
      for(int l = this->etdF(0).minL(); l < this->etdF(0).nOp(); ++l)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(l).constructBOperator(h*this->c(), basis.at(l).specLaplacian());

         // Compute the scaled exponential of the created operator
         this->computeScaledF0();

         // Compute the exponential from the possibly scaled exponential operator
         this->computeSquaredF0();

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rHarmOp(l).finaliseOp();
      }
   }

}
