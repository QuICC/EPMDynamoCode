/** \file ETD3Operators.cpp
 *  \brief Implementation of the ETD3 operators
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/Operators/ETD3Operators.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD3Operators::ETD3Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : ETDNOperators<4>(c, pTrunc, hasL0)
   {
   }

   void ETD3Operators::update(const EPMFloat dt, const ETD3Operators::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   void ETD3Operators::createOperators(const EPMFloat h, const ETD3Operators::BasisType &basis)
   {
      // Update the required scaling power
      this->updateScalings(h);

      // Loop over all degrees
      for(int l = this->etdF(0).minL(); l < this->etdF(0).nOp(); ++l)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(l).setOperator(h*this->c(), basis.at(l).specLaplacian().op());

         // Restrict the operator
         this->rEtdF(0).rHarmOp(l).constructBoundedOperator();

         // Compute the taylor expansion of the Fk functionals of the created operator
         this->computeScaledF3(l);

         // Compute the unscaled values
         this->computeSquaredF3(l);

         // Include the missing h factor
         this->rEtdF(1).rHarmOp(l).rOp() *= h;

         // Include the missing h factor
         this->rEtdF(2).rHarmOp(l).rOp() *= h;

         // Include the missing h^2 factor
         this->rEtdF(3).rHarmOp(l).rOp() *= h*h;

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rHarmOp(l).finaliseOp();
         this->rEtdF(1).rHarmOp(l).finaliseOp();
         this->rEtdF(2).rHarmOp(l).finaliseOp();
         this->rEtdF(3).rHarmOp(l).finaliseOp();
      }
   }

}
