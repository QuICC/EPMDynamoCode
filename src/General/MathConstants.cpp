/** \file MathConstants.cpp
 *  \brief Defines some useful math constants
 */

// System includes
//
#include <cmath>
#include <complex>

// External includes
//

// Class include
//
#include "General/MathConstants.hpp"

// Project includes
//

namespace EPMDynamo {

   const EPMFloat MathConstants::PI = std::acos(-1);

   const EPMComplex MathConstants::cI = EPMComplex(0.0, 1.0);

}
