/** \file InfluenceMatrix.cpp
 *  \brief Implementation of the influence matrix approach
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/InfluenceMatrix.hpp"

// Project includes
//

namespace EPMDynamo {

   InfluenceMatrix::InfluenceMatrix(SmartTruncation pTrunc, const InfluenceMatrix::BasisType &basis, bool hasL0)
      : LaplacianBOperatorSet<SimulationConfig::FactoredOpType>(basis, pTrunc, hasL0), mOpBCs(-2)
   {
      this->initSolutions();
   }

   void InfluenceMatrix::initSolutions()
   {
      int nN = this->trunc()->sim()->rad()->nN();
      int nL = this->trunc()->local()->spec()->nL();

      for(int l = 0; l < nL; ++l)
      {
         this->mSolutions.push_back(Array(nN));
         this->mSolutions.at(l).setConstant(0.0);
      }
   }

   void InfluenceMatrix::computeOperators()
   {
      // Create the laplacian operators
      this->createOperators(-1.0);
   }

}
