/** \file ETD2RKTraits.hpp
 *  \brief Traits to describe the requested ETD2RK method scheme
 */

#ifndef ETD2RKTRAITS_HPP
#define ETD2RKTRAITS_HPP

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
    * \brief Traits to some global aspects of the ETD2RK scheme
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD2RKTraits
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
         static const int INTERMEDIATE_STATES = 1;
   };

   template <typename TSimType> int ETD2RKTraits<TSimType>::order()
   {
      return 2;
   }; 

}

#endif // ETD2RKTRAITS_HPP
