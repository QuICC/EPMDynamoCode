/** \file TimestepControllerTraits.cpp
 *  \brief Soure of the timestep controller traits for the predefined controllers
 */

// System includes
//

// External includes
//

// Class includes
//
#include "Timestepping/Traits/TimestepControllerTraits.hpp"

// Project includes
//

namespace EPMDynamo {

   // 
   //  ELEMENTARY CONTROLLER
   //
   const EPMFloat TimestepControllerTraits<ElementaryCtrl>::KBETA1 = 1.0;

   const EPMFloat TimestepControllerTraits<ElementaryCtrl>::KBETA2 = 0.0;

   const EPMFloat TimestepControllerTraits<ElementaryCtrl>::ALPHA = 0.0;


   // 
   //  PI42 CONTROLLER
   //

   const EPMFloat TimestepControllerTraits<PI42Ctrl>::KBETA1 = 3.0/5.0;

   const EPMFloat TimestepControllerTraits<PI42Ctrl>::KBETA2 = -1.0/5.0;

   const EPMFloat TimestepControllerTraits<PI42Ctrl>::ALPHA = 0.0;


   // 
   //  H211B CONTROLLER
   //

   const EPMFloat TimestepControllerTraits<H211BCtrl>::KBETA1 = 1.0/4.0;

   const EPMFloat TimestepControllerTraits<H211BCtrl>::KBETA2 = 1.0/4.0;

   const EPMFloat TimestepControllerTraits<H211BCtrl>::ALPHA = 1.0/4.0;

}

