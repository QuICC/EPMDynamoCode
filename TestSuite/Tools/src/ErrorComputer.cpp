/** \file ErrorComputer.cpp
 *  This file contains the definitions of the ErrorComputer class methods
 */

// System includes
//

// External includes
//

// Class include
//
#include "ErrorComputer.hpp"

// Project includes
//

namespace EPMDynamo {

   ErrorComputer::ErrorComputer()
      : mMax(0.0), mRelMax(0.0)
   {
   }

   double ErrorComputer::RELATIVE_MAXERROR = 1.0e-10;

   double ErrorComputer::MAX_ALLOWED_ERROR = 1.0e-12;

   void ErrorComputer::resetErrors()
   {
      this->mMax = 0.0;

      this->mRelMax = 0.0;
   }

}
