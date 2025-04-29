/** \file ETD1Traits.hpp
 *  \brief Traits to describe the requested ETD1 method scheme
 */

#ifndef ETD1TRAITS_HPP
#define ETD1TRAITS_HPP

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
    * \brief Traits to some global aspects of the ETD1 scheme
    */
   class ETD1Traits
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
         static const int INTERMEDIATE_STATES = 0;
   };

}

#endif // ETD1TRAITS_HPP
