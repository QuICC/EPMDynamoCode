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
   template <typename TMethodTraits> class ETDTimestepControl: public TimestepControlBase
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
         /**
          * @brief Itermediate steps counter
          */
         int mInterCounter;
   };

   template <typename TMethodTraits> ETDTimestepControl<TMethodTraits>::ETDTimestepControl(TimestepParameters &tsParams, const EquationParameters &eqParams)
      : TimestepControlBase(tsParams, eqParams, TMethodTraits::CtrlType, TMethodTraits::order()), mInterCounter(0)
   {
   }

   template <typename TMethodTraits> void ETDTimestepControl<TMethodTraits>::checkConvergence(int step)
   {
      if(this->rTSParams().isNextStep())
      {
         // Switch to next steps in timestep iteration
         this->rTSParams().toggleTimestepStatus();

         // Accept timestep whatever happened before
         this->rTSParams().acceptTimestep();
      } else
      {
         // Check if timestep is bigger than minimum and if the initialisation steps are over
         if(step > TimestepConfig::MAX_INITIALISATION_STEPS && this->rTSParams().dt() < TimestepConfig::MIN_TIMESTEP)
         {
            // Trigger clean stop of simulation
            this->mKeepRunning = false;

            std::cout << "Timestep reached minimum timestep!" << std::endl;
         }
         // Check if Corrector correction norm is reducing
         else if(this->mInterCounter < TMethodTraits::INTERMEDIATE_STATES)
         {
            // Increment the intermediate steps counter
            ++this->mInterCounter;
         } else
         {
            // Reset number of intermediate steps counter
            this->mInterCounter = 0;

            // Toggle timestep to predictor step
            this->rTSParams().toggleTimestepStatus();

            // Reset timestep error
            this->resetError();
         }
      }
   }

   template <typename TMethodTraits> void ETDTimestepControl<TMethodTraits>::updateTimestep()
   {
      // Check timestep before starting new timestep iteration
      if(this->rTSParams().isNextStep())
      {
         // Store current timestep
         EPMFloat dt = this->rTSParams().dt();

         // Set timestep according to adaptive timestep controller
         if(! this->rTSParams().isRejected())
         {
            this->useAdaptiveTimestep(dt);
         }

         // Set timestep according to CFL condition
         this->testCFLCondition(dt);

         // Initialisation timestep
         this->testInitialisation(dt);

         // Force maximum timestep
         this->testMaximumTimestep(dt);

         // Use courant factor
         this->useCourantTimestep(dt);

         // Use windowed update scheme for timestep
         if(TMethodTraits::useWindowed)
         {
            this->useWindowedTimestep(dt);

         }
         // Use new timestep without window
         else
         {
            this->useTimestep(dt);
         }
      }
   }

}

#endif // ETDTIMESTEPCONTROL_HPP
