/** \file SimulationTruncation.cpp
 *  This file contains the definitions of the SimulationTruncation class methods
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/SimulationTruncation.hpp"

// Project includes
//

namespace EPMDynamo {

   SimulationTruncation::SimulationTruncation(const int maxN, const int nR, const int maxL, const int maxM, const int mp)
      : mpRadial(new RadialDomain(maxN, nR)), mpHorizon(new HorizontalDomain(maxL, maxM, mp))
   {
   }

}
