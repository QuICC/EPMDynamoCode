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
#include "Simulations/SimulationConstants.hpp"
#include "Config/SimulationConfig.hpp"
#include "Config/NumericalSchemeInc.hpp"

namespace EPMDynamo {

   RadialDomain::RadialDomain(const int maxN, const int nR)
      : RadialTruncation(maxN, nR), mMaxL(0), mpRSll(new Array(nR))
   {
   }

   void RadialDomain::setRadSll()
   {
      SimulationConfig::NumericalScheme::createRadSll(this->maxN(), this->mMaxL, this->radGrid(), this->mpRSll);
   }

}
