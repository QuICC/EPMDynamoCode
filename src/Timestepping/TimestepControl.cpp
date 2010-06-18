/** \file TimestepControl.cpp
 *  \brief Source of the implementation of a simple timestep control
 */

// System includes
//
#include <limits>
#include <iostream>

// Configuration includes
//
#include "Config/Parallelisation.h"

// External includes
//

// Class include
//
#include "Timestepping/TimestepControl.hpp"

// Project includes
//
#include "Timestepping/TimestepConfig.hpp"

namespace EPMDynamo {

   TimestepControl::TimestepControl(TimestepParameters &params, const EquationParameters &eqParams)
      : mKeepRunning(true), mNeedInit(true), mrParams(params), mrEqParams(eqParams), mCFLTimestep(0.0)
   {
   }

   void TimestepControl::updateCFLTimestep(const RTPField &velV, const EPMFloat llFactor)
   {
      if(this->mrParams.isNextStep())
      {
         // Get truncation information
         int r0 = velV.trunc()->local()->rtp()->r0();
         int nR = velV.nR();
         int nRSim = velV.trunc()->sim()->rad()->nR();

         // Get radii
         Array radii = velV.trunc()->sim()->rad()->radGrid();

         // Initialise to highest representable value
         this->mCFLTimestep = std::numeric_limits<EPMFloat>::max();

         // Local CFL Condition
         EPMFloat dr;
         EPMFloat maxVel = 0.0;
         Matrix   p;
         int n_;
         for(int n=0; n < nR; ++n)
         {
            n_ = n + r0;
            if(n_ == 0)
            {
               dr = radii(1)-radii(0);
            }
            else if(n_ == nRSim - 1)
            {
               dr = 1.0 - radii(n_);
            } else
            {
               dr = std::min(radii(n_) - radii(n_-1), radii(n_+1) - radii(n_));
            }

            // Radial Velocity
            maxVel = velV.r().shell(n).cwise().abs().maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               this->mCFLTimestep = std::min(this->mCFLTimestep, dr/maxVel);
            }

            // Angular "Velocity"
            dr = radii(n_)/std::sqrt(llFactor);
            maxVel = (velV.theta().shell(n).cwise().square() + velV.phi().shell(n).cwise().square()).cwise().sqrt().maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               this->mCFLTimestep = std::min(this->mCFLTimestep, dr/maxVel);
            }
         }

         // Global CFL conditions min(Ro, sqrt(E))
         if(this->mrEqParams.Ro() != 0.0)
         {
            this->mCFLTimestep = std::min(this->mCFLTimestep, this->mrEqParams.Ro());
         }

         if(this->mrEqParams.E() != 0.0)
         {
            this->mCFLTimestep = std::min(this->mCFLTimestep, std::sqrt(this->mrEqParams.E()));
         }
      }
   }

   void TimestepControl::updateCFLTimestep(const RTPField &magB, const RTPField &velV, const EPMFloat llFactor)
   {
      if(mrParams.isNextStep())
      {
         // Get truncation information
         int r0 = magB.trunc()->local()->rtp()->r0();
         int nR = magB.nR();
         int nRSim = magB.trunc()->sim()->rad()->nR();

         // Get radii
         Array radii = magB.trunc()->sim()->rad()->radGrid();

         // Initialise to highest representable value
         this->mCFLTimestep = std::numeric_limits<EPMFloat>::max();

         // Local CFL Condition
         EPMFloat d;
         EPMFloat dr;
         EPMFloat maxVel = 0.0;
         Matrix   p;
         int n_;
         for(int n=0; n < nR; ++n)
         {
            n_ = n + r0;
            if(n_ == 0)
            {
               dr = radii(1)-radii(0);
            }
            else if(n_ == nRSim - 1)
            {
               dr = 1.0 - radii(n_);
            } else
            {
               dr = std::min(radii(n_) - radii(n_-1), radii(n_+1) - radii(n_));
            }
            d = (mrEqParams.E() + mrEqParams.Ro())/(2.0*dr);
            d = d*d;

            // Radial "Velocity"
            p = magB.r().shell(n).cwise().square();
            maxVel = (p.cwise()/((p*mrEqParams.Ro()).cwise() + d).cwise().sqrt() + velV.r().shell(n).cwise().abs()).maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               this->mCFLTimestep = std::min(this->mCFLTimestep, dr/maxVel);
            }

            // Angular "Velocity"
            dr = radii(n_)/std::sqrt(llFactor);
            d = (mrEqParams.E() + mrEqParams.Ro())/(2.0*dr);
            d = d*d;
            p = magB.theta().shell(n).cwise().square() + magB.phi().shell(n).cwise().square();
            maxVel = (p.cwise()/((p*mrEqParams.Ro()).cwise() + d).cwise().sqrt() + (velV.theta().shell(n).cwise().square() + velV.phi().shell(n).cwise().square()).cwise().sqrt()).maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               this->mCFLTimestep = std::min(this->mCFLTimestep, dr/maxVel);
            }
         }

         // Global CFL conditions min(Ro, sqrt(E))
         if(this->mrEqParams.Ro() != 0.0)
         {
            this->mCFLTimestep = std::min(this->mCFLTimestep, this->mrEqParams.Ro());
         }

         if(this->mrEqParams.E() != 0.0)
         {
            this->mCFLTimestep = std::min(this->mCFLTimestep, std::sqrt(this->mrEqParams.E()));
         }
      }
   }

   void TimestepControl::getSimulationCFLCondition()
   {
      // Get the "global" local minimum for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mCFLTimestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   void TimestepControl::checkConvergence(int step)
   {
      if(this->mrParams.isNextStep())
      {
         // Switch step type
         this->mrParams.toggleTimestepStatus();
      } else
      {
         // Check if timestep value after initialisation steps
         if(this->mrParams.dt() < TimestepConfig::MIN_TIMESTEP)
         {
            // Trigger clean stop of simulation
            this->mKeepRunning = false;

            std::cout << "Timestep reached minimum timestep!" << std::endl;
         }
         // Timestepping converged next timestep will be a predictor one
         else
         {
            // Switch step type
            this->mrParams.toggleTimestepStatus();
         }
      }
   }

   void TimestepControl::testInitialisationTimestep(EPMFloat& rDt)
   {
      // Initialisation timestep
      if(this->mNeedInit)
      {
         // Set simple starting timestep to mininal timestep
         if(rDt == TimestepConfig::MIN_TIMESTEP)
         {
            rDt = TimestepConfig::MIN_TIMESTEP*2.0;
         } else
         {
            rDt = std::min(rDt/1.0e3, TimestepConfig::MIN_TIMESTEP*2.0);
         }

         std::cout << "Use Initialisation Timestep!" << std::endl;

         // Set flag to forbid second entry
         this->mNeedInit = false;
      }
   }

   void TimestepControl::testMaximumTimestep(EPMFloat& rDt)
   {
      // Force maximum timestep
      if(rDt > TimestepConfig::MAXIMUM_TIMESTEP)
      {
         // Set timestep to maximum allowed timestep
         rDt = TimestepConfig::MAXIMUM_TIMESTEP;
      }
   }

   void TimestepControl::testCFLTimestep(EPMFloat& rDt)
   {
      // Set timestep according to CFL condition
      this->getSimulationCFLCondition();
      if(this->mCFLTimestep > 0.0)
      {
         rDt = std::min(this->mCFLTimestep, rDt);
      }
   }

   void TimestepControl::useCourantTimestep(EPMFloat& rDt)
   {
      // Include Courant number into timestep value
      rDt = rDt*TimestepConfig::COURANT_NUMBER;
   }

   void TimestepControl::setWindowedTimestep(EPMFloat dt)
   {
      // Don't make adaptative scheme to reactive. Allow a range of timestep values
      if(dt < this->mrParams.dt()*(1.0-TimestepConfig::TIMESTEP_MARGIN) ||  dt > this->mrParams.dt()*(1.0 + TimestepConfig::TIMESTEP_MARGIN))
      {
         // Store new timestep
         this->mrParams.updateTimestep(dt);

         std::cout << "Updating Timestep! Now using dt=" << this->mrParams.dt() << " at time t=" << this->mrParams.time()<< std::endl;
      }
   }

   void TimestepControl::updateTimestep()
   {
      if(this->mrParams.isNextStep())
      {
         // Get current dt
         EPMFloat dt = this->mrParams.dt();

         // Set timestep according to CFL condition
         this->testCFLTimestep(dt);

         // Initialisation timestep
         this->testInitialisationTimestep(dt);

         // Force maximum timestep
         this->testMaximumTimestep(dt);

         // Use courant factor
         this->useCourantTimestep(dt);

         // Don't make adaptative scheme to reactive. Allow a range of timestep values
         this->setWindowedTimestep(dt);
      }
   }

}
