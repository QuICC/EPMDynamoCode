/** \file EquationParameters.cpp
 *  This file contains the definitions of the EquationParameters class methods
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

   EquationParameters::EquationParameters(DynamoFloat E, DynamoFloat q, DynamoFloat Ra, DynamoFloat Ro)
      : mE(E), mQ(q), mRa(Ra), mRo(Ro)
   {
   }

   EquationParameters::EquationParameters(const Array& arr)
      : mE(arr(0)), mQ(arr(1)), mRa(arr(2)), mRo(arr(3))
   {
   }

}
