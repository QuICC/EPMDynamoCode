/** \file ThetaTraits.cpp
 *  \brief Traits to describe the requested PC \f$\theta\f$-method scheme
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/PredictorCorrector/Theta/ThetaTraits.hpp"

// Project includes
//

namespace EPMDynamo {

   const EPMFloat ThetaTraits::theta = 0.5; 

   int ThetaTraits::order()
   {
      if(ThetaTraits::theta == 0.5)
      {
         return 2;
      } else
      {
         return 1;
      }
   }; 

}
