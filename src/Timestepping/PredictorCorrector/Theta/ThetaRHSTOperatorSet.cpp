/** \file ThetaRHSTOperatorSet.cpp
 *  \brief Implementation of the \f$\theta\f$ method RHS operator set
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class includes
//
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"

// Project includes
//

namespace EPMDynamo {

   ThetaRHSTOperatorSet::ThetaRHSTOperatorSet(EPMFloat a, EPMFloat b, const ThetaRHSTOperatorSet::BasisType &basis, SmartTruncation pTrunc, bool hasL0)
      : ThetaTOperatorSet<SimulationConfig::OperatorType, ExplicitTOperatorSet>(a, b, basis, pTrunc, hasL0)
   {
   }

}
