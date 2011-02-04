/** \file RadialDomain.cpp
 *  \brief Source of the implementation of the radial domain information
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Simulation/RadialDomain.hpp"

// Project includes
//

namespace EPMDynamo {

   RadialDomain::RadialDomain(const int maxN, const int nR)
      : RadialTruncation(maxN, nR), mpRSll(new Array(nR))
   {
   }

}
