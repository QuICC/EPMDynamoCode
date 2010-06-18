/** \file ETDTimestepControl.cpp
 *  This file contains the definitions of the ETDTimestepControl class methods
 */

// System includes
//
#include <limits>
#include <iostream>
#include "Config/Parallelisation.h"

// External includes
//

// Class include
//
#include "Timestepping/ETD/ETDTimestepControl.hpp"

// Project includes
//

namespace EPMDynamo {

   ETDTimestepControl::ETDTimestepControl(TimestepParameters &params, const EquationParameters &eqParams)
      : TimestepControl(params, eqParams)
   {
   }

   ETDTimestepControl::~ETDTimestepControl()
   {
   }

}
