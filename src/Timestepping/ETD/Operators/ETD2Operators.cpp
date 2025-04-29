/** \file ETD2Operators.cpp
 *  \brief Implementation of the ETD2 operators
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/Operators/ETD2Operators.hpp"

// Project includes
//

namespace EPMDynamo {

   ETD2Operators::ETD2Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : ETDNOperators<3>(c, pTrunc, hasL0)
   {
   }

   void ETD2Operators::update(const EPMFloat dt, const ETD2Operators::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   void ETD2Operators::createOperators(const EPMFloat h, const ETD2Operators::BasisType &basis)
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
         this->computeScaledF2(l);

         // Compute the unscaled values
         this->computeSquaredF2(l);

         // Include the missing h factor
         this->rEtdF(1).rHarmOp(l).rOp() *= h;

         // Include the missing h factor
         this->rEtdF(2).rHarmOp(l).rOp() *= h;

         // Correct the results for zero eigenvalues
         if(!this->isFullRank(l))
         { 
            this->rEtdF(1).rHarmOp(l).rOp().topRows(1) += h*this->etdF(0).harmOp(l).op().leftCols(1).transpose();

            this->rEtdF(2).rHarmOp(l).rOp().topRows(1) += 0.5*h*this->etdF(0).harmOp(l).op().leftCols(1).transpose();
         }

         // Do finalisation step
         this->rEtdF(0).rHarmOp(l).finaliseOp();
         this->rEtdF(1).rHarmOp(l).finaliseOp();
         this->rEtdF(2).rHarmOp(l).finaliseOp();
      }
   }

}
