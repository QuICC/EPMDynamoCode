/** \file ThetaTraits.hpp
 *  \brief Traits to describe the requested PC \f$\theta\f$-method scheme
 */

#ifndef THETATRAITS_HPP
#define THETATRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   // Forward declarations
   template <typename> class ErrorL2Max;

   /**
    * \brief Traits to some global aspects of the theta scheme
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaTraits
   {
      public:
         /// Implicitness value of the theta scheme
         static const EPMFloat theta;

         /// Typedef for the error norm type
         typedef ErrorL2Max<TSimType>   ErrorNormType;

         /// Type of timestep controller to use
         static const TimestepCtrlTypes  CtrlType = PI42Ctrl;

         /// Order of the timestep scheme
         static int order();

         /// Use windowes adaptive timestep
         static const bool useWindowed = false;
   };

   template <typename TSimType> const EPMFloat ThetaTraits<TSimType>::theta = 0.5; 

   template <typename TSimType> int ThetaTraits<TSimType>::order()
   {
      if(ThetaTraits<TSimType>::theta == 0.5)
      {
         return 2;
      } else
      {
         return 1;
      }
   }; 

}

#endif // THETATRAITS_HPP
