/** \file ETDTimestepControl.hpp
 *  \brief Implementation of a specialised timestep control for ETD schemes
 */

#ifndef ETDTIMESTEPCONTROL_HPP
#define ETDTIMESTEPCONTROL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "Timestepping/TimestepControlBase.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/TimestepConfig.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a specialised timestep control for predictor/corrector schemes
    */
   template <typename TMethodTraits> class ETDTimestepControl: public TimestepControl
   {
      public:
         /**
          * @brief Constructor
          *
          * @param tsParams Timestep parameters
          * @param eqParams Equation parameters
          */
         ETDTimestepControl(TimestepParameters &tsParams, const EquationParameters &eqParams);

         /**
          * @brief Destructor
          */
         virtual ~ETDTimestepControl() {};

         /**
          * @brief Check that the timestepping is converging
          *
          * @param step Current step
          */
         virtual void checkConvergence(int step);

         /**
          * @brief Update the timestep to use
          */
         virtual void updateTimestep();
         
      protected:

      private:
   };

   template <typename TMethodTraits> ETDTimestepControl<TMethodTraits>::ETDTimestepControl(TimestepParameters &tsParams, const EquationParameters &eqParams)
      : TimestepControlBase(tsParams, eqParams, TMethodTraits::CtrlType, TMethodTraits::order())
   {
   }

   template <typename TMethodTraits> void ETDTimestepControl<TMethodTraits>::checkConvergence(int step)
   {
   }

   template <typename TMethodTraits> void ETDTimestepControl<TMethodTraits>::updateTimestep()
   {
   }

}

#endif // ETDTIMESTEPCONTROL_HPP
