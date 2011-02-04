/** \file SimulationTruncation.cpp
 *  \brief Source of the Simulation wide truncation information
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
      // Initialise the radial sqrt(L(L+1)) factors
      this->mpRadial->initRadSll(maxL);

      // FOR THE MOMENT THE VALUES ARE SET BY HAND
      this->mpRadial->setRadSll();
   }

}
