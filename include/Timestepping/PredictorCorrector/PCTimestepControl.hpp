/** \file PCTimestepControl.hpp
 *  \brief Implementation of a specialised timestep control for predictor/corrector schemes
 */

#ifndef PCTIMESTEPCONTROL_HPP
#define PCTIMESTEPCONTROL_HPP

// System includes
//
#include <limits>
#include <iostream>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "Timestepping/TimestepControlBase.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/TimestepConfig.hpp"
#include "Timestepping/PredictorCorrector/PCTimestepConfig.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a specialised timestep control for predictor/corrector schemes
    */
   template <typename TMethodTraits> class PCTimestepControl: public TimestepControlBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param tsParams Timestep parameters
          * @param eqParams Equation parameters
          */
         PCTimestepControl(TimestepParameters &tsParams, const EquationParameters &eqParams);

         /**
          * @brief Constructor
          */
         virtual ~PCTimestepControl() {};

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
         /**
          * @brief Counter for the number of Corrector iterations
          */
         int  mCorrectorIts;

         /**
          * @brief Storage for the norm of the current corrector iteration
          */
         EPMFloat   mCorrectionNorm;

         /**
          * @brief Storage for the norm of the previous corrector iteration
          */
         EPMFloat   mOldCorrectionNorm;

         /**
          * @brief Get PC errors and compute required timestep
          */
         void setErrorConstrains();

      private:
   };

   template <typename TMethodTraits> PCTimestepControl<TMethodTraits>::PCTimestepControl(TimestepParameters &tsParams, const EquationParameters &eqParams)
      : TimestepControlBase(tsParams, eqParams, TMethodTraits::CtrlType, TMethodTraits::order()), mCorrectorIts(0), mCorrectionNorm(0.0), mOldCorrectionNorm(0.0)
   {
   }

   template <typename TMethodTraits> void PCTimestepControl<TMethodTraits>::setErrorConstrains()
   {
      // Get Corrector correction norm
      this->mCorrectionNorm = this->rTSParams().error();

      // Reset Corrector correction norm
      this->rTSParams().resetError();

      // Use total corrector correction
      this->mError += this->mCorrectionNorm;
   }

   template <typename TMethodTraits> void PCTimestepControl<TMethodTraits>::checkConvergence(int step)
   {
      if(this->rTSParams().isNextStep())
      {
         // Switch to corrector steps
         this->rTSParams().toggleTimestepStatus();

         // Prepare for first corrector step 
         this->mCorrectorIts = 1;

         // Accept timestep whatever happened before
         this->rTSParams().acceptTimestep();
      } else
      {
         // Compute the error norm timestep and reset error
         this->setErrorConstrains();

         // Initialise previous correction norm to maximum value
         if(this->mCorrectorIts == 1)
         {
            // Initialise previous corrector norm to max representable value
            this->mOldCorrectionNorm = std::numeric_limits<EPMFloat>::max();
         }

         // Check if timestep is bigger than minimum and if the initialisation steps are over
         if(step > TimestepConfig::MAX_INITIALISATION_STEPS && this->rTSParams().dt() < TimestepConfig::MIN_TIMESTEP)
         {
            // Reset iteration counter
            this->mCorrectorIts = 0;

            // Trigger clean stop of simulation
            this->mKeepRunning = false;

            std::cout << "Timestep reached minimum timestep!" << std::endl;
         }
         // Timestepping didn't converge and reached max number of allowed corrector iterations
         else if(this->mCorrectorIts > PCTimestepConfig::MAX_CORRECTOR_ITERATIONS)
         {
            // Reset iteration counter
            this->mCorrectorIts = 0;

            // Trigger clean stop of simulation
            this->mKeepRunning = false;

            std::cout << "Timestep reached max corrector iterations!" << std::endl;
         }
         // Check if Corrector correction norm is reducing
         else if(this->mCorrectionNorm > this->mOldCorrectionNorm)       
         {
            // Stop computation if corrector norm is more than twice the max value
            if(this->mCorrectionNorm > 2.0*PCTimestepConfig::MAX_CORRECTION_NORM)
            {
               // Trigger clean stop of simulation
               this->mKeepRunning = false;

               std::cout << "Corrector norm is increasing!" << std::endl;
            // .. else give it a last try with a 10x smaller timestep
            } else
            {
               // Print information on output
               std::cout << "Timestep control failed, rejecting timestep!" << std::endl;

               // Try again with 10x smaller timestep
               this->useTimestep(this->rTSParams().dt()/(10.0*TimestepConfig::COURANT_NUMBER));

               // Trigger reject of previous timestep (reinitialises fields to previous values)
               this->rTSParams().rejectTimestep();

               // Toggle timestep to predictor step
               this->rTSParams().toggleTimestepStatus();
            }

            // Reset corrector iterations counter
            this->mCorrectorIts = 0;
         }
         // Check if corrector error is bigger than threshold
         else if(this->mError > PCTimestepConfig::MAX_CORRECTOR_ERROR)
         {
            // Reset iteration counter
            this->mCorrectorIts = 0;

            // Trigger clean stop of simulation
            this->mKeepRunning = false;

            std::cout << "Timestep reached max corrector error!" << std::endl;
         }
         // Check if corrector correction norm is bigger than threshold
         else if(this->mCorrectionNorm > PCTimestepConfig::MAX_CORRECTION_NORM)
         {
            // Update old corrector correction norm value
            this->mOldCorrectionNorm = this->mCorrectionNorm;

            // Increment corrector step counter
            ++this->mCorrectorIts;

            // Output some information
            if(this->mCorrectorIts == 2)
            {
               std::cout << "---> Requires multiple correctors:" << std::endl;
            }
            std::cout << "------> Corrector: " << this->mCorrectorIts << std::endl;

         }
         // Timestepping converged next timestep will be a predictor one
         else
         {
            // reset corrector iterations counter
            this->mCorrectorIts = 0;

            // Toggle timestep to predictor step
            this->rTSParams().toggleTimestepStatus();
         }
      }
   }

   template <typename TMethodTraits> void PCTimestepControl<TMethodTraits>::updateTimestep()
   {
      // Check if predictor step is coming
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

         // Use new timestep without window
         } else
         {
            this->useTimestep(dt);
         }

         // Reset timestep error
         this->resetError();
      }
   }

}

#endif // PCTIMESTEPCONTROL_HPP
