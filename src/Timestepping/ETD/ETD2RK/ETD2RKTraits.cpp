/** \file ETD2RKTraits.cpp
 *  \brief Traits to describe the requested ETD2RK method scheme
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/ETD2RK/ETD2RKTraits.hpp"

// Project includes
//

namespace EPMDynamo {

   int ETD2RKTraits::order()
   {
      return 2;
   }; 

}
