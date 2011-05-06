/** \file ETD2RKTraits.hpp
 *  \brief Traits to describe the requested ETD2RK method scheme
 */

#ifndef ETD2RKTRAITS_HPP
#define ETD2RKTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits to some global aspects of the ETD2RK scheme
    */
   class ETD2RKTraits
   {
      public:
         /// Typedef for the error norm type
         typedef ErrorModeL2Max   ErrorNormType;

         /// Type of timestep controller to use
         static const TimestepCtrlTypes  CtrlType = H211BCtrl;

         /// Order of the timestep scheme
         static int order();

         /// Use windowes adaptive timestep
         static const bool useWindowed = false;

         /// Number of intermediate states
         static const int INTERMEDIATE_STATES = 1;
   };

}

#endif // ETD2RKTRAITS_HPP
