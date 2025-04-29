/** \file ErrorControl.cpp
 *  \brief Implementation of a specialised error control
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ErrorControl.hpp"

// Project includes
//

namespace EPMDynamo {

   EPMFloat ErrorControl::errorNorm(const ErrorControl::ScalarType& rVar,const ErrorControl::ScalarType& rRef, const EPMFloat oldError)
   {
      // Create temporary storage
      EPMFloat norm = 0.0;
      EPMFloat error = 0.0;

      // Compute error norm
      norm = ErrorNormType::computeNorm(rVar, rRef);

      // Update error norm
      error = ErrorNormType::updateNorm(norm, oldError);

      return error;
   }

}
