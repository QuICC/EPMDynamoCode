/** \file ETD1Operators.cpp
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
#include "Timestepping/ETD/Operators/ETD1Operators.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD1Operators::ETD1Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : ETDNOperators<2>(c, pTrunc, hasL0)
   {
   }

   void ETD1Operators::update(const EPMFloat dt, const ETD1Operators::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   void ETD1Operators::createOperators(const EPMFloat h, const ETD1Operators::BasisType &basis)
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
         this->computeScaledF1(l);

         // Compute the unscaled values
         this->computeSquaredF1(l);

         // Include the missing h factor
         this->rEtdF(1).rHarmOp(l).rOp() *= h;

         // Correct the results for zero eigenvalues
         if(!this->isFullRank(l))
         { 
            this->rEtdF(1).rHarmOp(l).rOp().topRows(1) += h*this->etdF(0).harmOp(l).op().leftCols(1).transpose();
         }

         // Do finalisation step
         this->rEtdF(0).rHarmOp(l).finaliseOp();
         this->rEtdF(1).rHarmOp(l).finaliseOp();
      }
   }

}
