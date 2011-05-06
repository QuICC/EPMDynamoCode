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
#include "Timestepping/ErrorNorms/ErrorModeL2Max.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits to some global aspects of the theta scheme
    */
   class ThetaTraits
   {
      public:
         /// Implicitness value of the theta scheme
         static const EPMFloat theta;

         /// Typedef for the error norm type
         typedef ErrorModeL2Max   ErrorNormType;

         /// Type of timestep controller to use
         static const TimestepCtrlTypes  CtrlType = H211BCtrl;

         /// Order of the timestep scheme
         static int order();

         /// Use windowes adaptive timestep
         static const bool useWindowed = false;
   };

}

#endif // THETATRAITS_HPP
