/** \file ThetaLHSTOperatorSet.cpp
 *  \brief Implementation of the \f$\theta\f$-method LHS operator set
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"

// Project includes
//

namespace EPMDynamo {

   ThetaLHSTOperatorSet::ThetaLHSTOperatorSet(EPMFloat a, EPMFloat b, const ThetaLHSTOperatorSet::BasisType &basis, SmartTruncation pTrunc, bool hasL0)
      : ThetaTOperatorSet<SimulationConfig::FactoredOpType, ImplicitTOperatorSet>(a, b, basis, pTrunc, hasL0)
   {
   }

}
