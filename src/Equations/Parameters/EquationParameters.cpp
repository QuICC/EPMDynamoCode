/** \file EquationParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/EquationParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   EquationParameters::EquationParameters(EPMFloat E, EPMFloat q, EPMFloat Ra, EPMFloat Ro)
      : mE(E), mQ(q), mRa(Ra), mRo(Ro)
   {
   }

   EquationParameters::EquationParameters(const Array& arr)
      : mE(arr(0)), mQ(arr(1)), mRa(arr(2)), mRo(arr(3))
   {
   }

}
