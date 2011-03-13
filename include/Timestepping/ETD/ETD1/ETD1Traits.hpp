/** \file ETD1Traits.hpp
 *  \brief Traits to describe the requested ETD1 method scheme
 */

#ifndef ETD1TRAITS_HPP
#define ETD1TRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   // Forward declarations
   template <typename> class ErrorModeL2Max;

   /**
    * \brief Traits to some global aspects of the ETD1 scheme
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD1Traits
   {
      public:
         /// Typedef for the error norm type
         typedef ErrorModeL2Max<TSimType>   ErrorNormType;

         /// Type of timestep controller to use
         static const TimestepCtrlTypes  CtrlType = H211BCtrl;

         /// Order of the timestep scheme
         static int order();

         /// Use windowes adaptive timestep
         static const bool useWindowed = false;

         /// Number of intermediate states
         static const int INTERMEDIATE_STATES = 0;
   };

   template <typename TSimType> int ETD1Traits<TSimType>::order()
   {
      return 1;
   }; 

}

#endif // ETD1TRAITS_HPP
