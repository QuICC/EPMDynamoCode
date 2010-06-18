/** \file PCTimestepControl.cpp
 *  \brief Source of the implementation of the specialised timestep control for predictor/corrector schemes
 */

// System includes
//
#include <limits>
#include <iostream>

// External includes
//

// Class include
//
#include "Timestepping/PredictorCorrector/PCTimestepControl.hpp"

// Project includes
//
#include "Timestepping/TimestepConfig.hpp"
#include "Timestepping/PredictorCorrector/PCTimestepConfig.hpp"

namespace EPMDynamo {

   PCTimestepControl::PCTimestepControl(TimestepParameters &params, const EquationParameters &eqParams)
      : TimestepControl(params, eqParams), mCounterCorrIts(0), mPCTStep(0.0), mCorrectionNorm(0.0), mOldCorrectionNorm(0.0), mPCError(TimestepConfig::TIMESTEP_ERROR_EPSILON), mOldPCError(TimestepConfig::TIMESTEP_ERROR_EPSILON), mController(TimestepController::ELEMENTARY, 2, params), mTotCorr(0.0)
   {
   }

   void PCTimestepControl::setErrorConstrains()
   {
      // Get Corrector correction norm
      this->mCorrectionNorm = this->mrParams.error();

      // Reset Corrector correction norm
//      this->mrParams.resetError();

this->mTotCorr += this->mCorrectionNorm;

      // Update PC error diagnostic value
      if(PCTimestepConfig::PC_ERROR_TYPE == PCTimestepConfig::CORRECTOR_NORM)
      {
         // Use last corrector correction value
         this->mPCError = this->mCorrectionNorm;

      } else if(PCTimestepConfig::PC_ERROR_TYPE == PCTimestepConfig::SUMMED_CORRECTOR_NORM)
      {
         // Use total corrector correction
         this->mPCError += this->mCorrectionNorm;
      }
   }

   void PCTimestepControl::computePCTimestep()
   {
      this->mPCTStep = this->mController.nextTimestep(this->mPCError, this->mOldPCError);
   }

   void PCTimestepControl::resetPCError()
   {
      // Store PC error as old value
      this->mOldPCError = this->mPCError;
      
      // Reset PC error to zero
      this->mPCError = 0.0;

this->mTotCorr = 0.0;
   }

   void PCTimestepControl::checkConvergence(int step)
   {
      if(this->mrParams.isNextStep())
      {
         // Switch to corrector steps
         this->mrParams.toggleTimestepStatus();

         // Prepare for first corrector step 
         this->mCounterCorrIts = 1;

         // Accept timestep whatever happened before
         this->mrParams.acceptTimestep();
      } else
      {
         // Compute the error norm timestep and reset error
         this->setErrorConstrains();

         // Initialise previous correction norm to maximum value
         if(this->mCounterCorrIts == 1)
         {
            // Initialise previous corrector norm to max representable value
            this->mOldCorrectionNorm = std::numeric_limits<EPMFloat>::max();
         }

         // Check if timestep is bigger than minimal value if the initialisation steps are over
         if(step > PCTimestepConfig::MAX_INITIALISATION_STEPS && this->mrParams.dt() < TimestepConfig::MIN_TIMESTEP)
         {
            // Reset iteration counter
            this->mCounterCorrIts = 0;

            // Trigger clean stop of simulation
            this->mKeepRunning = false;

            std::cout << "Timestep reached minimum timestep!" << std::endl;
         }
         // Timestepping didn't converged and reached max number of iterations
         else if(this->mCounterCorrIts > PCTimestepConfig::MAX_CORRECTOR_ITERATIONS)
         {
            // Reset iteration counter
            this->mCounterCorrIts = 0;

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
               // Try again with 10x smaller timestep
               this->mPCTStep = this->mrParams.dt()/(10.0*TimestepConfig::COURANT_NUMBER);

               // Trigger reject of previous timestep (reinitialises fields to previous values)
               this->mrParams.rejectTimestep();

               // Toggle timestep to predictor step
               this->mrParams.toggleTimestepStatus();

               // Print information on output
               std::cout << "Timestep control failed, rejecting timestep!" << std::endl;
            }

            // Reset corrector iterations counter
            this->mCounterCorrIts = 0;
         }
         // Check if corrector correction norm is bigger than threshold
         else if(this->mCorrectionNorm > PCTimestepConfig::MAX_CORRECTION_NORM)
         {
            // Update old corrector correction norm value
            this->mOldCorrectionNorm = this->mCorrectionNorm;

            // Increment corrector step counter
            ++this->mCounterCorrIts;
         }
         // Timestepping converged next timestep will be a predictor one
         else
         {
            std::cout << "Converged after " << this->mCounterCorrIts << " Tot corr: " << this->mTotCorr << " err: " << this->mPCError << std::endl;

            // reset corrector iterations counter
            this->mCounterCorrIts = 0;

            // Compute corrector timestep
            this->computePCTimestep();

            // Toggle timestep to predictor step
            this->mrParams.toggleTimestepStatus();

            // Reset corrector norm
            this->resetPCError();
         }
      }
   }

   void PCTimestepControl::testPCTimestep(EPMFloat& rDt)
   {
      // Set timestep according to PC errors
      if(this->mPCTStep > 0.0)
      {
         // Forbid big timestep increases (limit to a factor 2)
         //rDt = std::min(this->mPCTStep, 2.0*rDt/TimestepConfig::COURANT_NUMBER);
         
         // Give full control to TimestepController
         rDt = this->mPCTStep;
      }
   }

   void PCTimestepControl::updateTimestep()
   {
      // Check if predictor step is coming
      if(mrParams.isNextStep())
      {
         // Store current timestep
         EPMFloat dt = this->mrParams.dt();

         // Set timestep according to corrector correction norm
         this->testPCTimestep(dt);

         // Set timestep according to CFL condition
         this->testCFLTimestep(dt);

         // Initialisation timestep
         this->testInitialisationTimestep(dt);

         // Force maximum timestep
         this->testMaximumTimestep(dt);

         // Use courant factor
         this->useCourantTimestep(dt);

         // Don't make adaptative scheme to reactive. Allow a range of timestep values
         //this->setWindowedTimestep(dt);
         
         this->mrParams.updateTimestep(dt);
      }
   }

}
