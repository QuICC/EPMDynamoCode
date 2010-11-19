/** \file TimestepControlBase.cpp
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
#include "Timestepping/TimestepControlBase.hpp"

// Project includes
//
#include "Timestepping/TimestepConfig.hpp"

namespace EPMDynamo {

   TimestepControlBase::TimestepControlBase(TimestepParameters &tsParams, const EquationParameters &eqParams, TimestepCtrlTypes ctrlType, int order)
      : mKeepRunning(true), mNeedInit(true), mError(TimestepConfig::TIMESTEP_ERROR_EPSILON), mOldError(TimestepConfig::TIMESTEP_ERROR_EPSILON), mCFLTimestep(0.0), mrTSParams(tsParams), mrEqParams(eqParams), mController(ctrlType, order, tsParams)
   {
   }

   void TimestepControlBase::resetError()
   {
      // Store error as old value
      this->mOldError = this->mError;
      
      // Reset error to zero
      this->mError = 0.0;
   }

   void TimestepControlBase::updateCFLTimestep(const RTPField &velV, const EPMFloat llFactor)
   {
      if(this->rTSParams().isNextStep())
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
            maxVel = velV.r().shell(n).array().abs().maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               this->mCFLTimestep = std::min(this->mCFLTimestep, dr/maxVel);
            }

            // Angular "Velocity"
            dr = radii(n_)/std::sqrt(llFactor);
            maxVel = (velV.theta().shell(n).array().square() + velV.phi().shell(n).array().square()).array().sqrt().maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               this->mCFLTimestep = std::min(this->mCFLTimestep, dr/maxVel);
            }
         }

         // Global CFL conditions min(Ro, sqrt(E))
         this->eqParams().testGlobalCFL(this->mCFLTimestep);
      }
   }

   void TimestepControlBase::updateCFLTimestep(const RTPField &magB, const RTPField &velV, const EPMFloat llFactor)
   {
      if(this->rTSParams().isNextStep())
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
            d = this->eqParams().alfvenDamping(dr);

            // Radial "Velocity"
            p = magB.r().shell(n).array().square()*this->eqParams().alfvenFactor();
            maxVel = (p.array()/(p.array() + d).array().sqrt() + velV.r().shell(n).array().abs()).maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               this->mCFLTimestep = std::min(this->mCFLTimestep, dr/maxVel);
            }

            // Angular "Velocity"
            dr = radii(n_)/std::sqrt(llFactor);
            d = this->eqParams().alfvenDamping(dr);
            p = (magB.theta().shell(n).array().square() + magB.phi().shell(n).array().square())*this->eqParams().alfvenFactor();
            maxVel = (p.array()/(p.array() + d).array().sqrt() + (velV.theta().shell(n).array().square() + velV.phi().shell(n).array().square()).array().sqrt()).maxCoeff();

            // Update timestep
            if(maxVel != 0.0)
            {
               this->mCFLTimestep = std::min(this->mCFLTimestep, dr/maxVel);
            }
         }

         // Global CFL conditions min(Ro, sqrt(E))
         this->eqParams().testGlobalCFL(this->mCFLTimestep);
      }
   }

   void TimestepControlBase::getSimulationCFLCondition()
   {
      // Get the "global" local minimum for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mCFLTimestep, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   void TimestepControlBase::testInitialisation(EPMFloat& rDt)
   {
      // Initialisation timestep
      if(this->mNeedInit)
      {
         // Set simple starting timestep to mininal timestep
         if(rDt == TimestepConfig::MIN_TIMESTEP)
         {
            rDt = TimestepConfig::MIN_TIMESTEP*TimestepConfig::TIMESTEP_MAX_RATIO;
         } else
         {
            rDt = std::max(rDt/1.0e3, TimestepConfig::MIN_TIMESTEP*TimestepConfig::TIMESTEP_MAX_RATIO);
         }

         // Set flag to forbid second entry
         this->mNeedInit = false;

         std::cout << "Use Initialisation Timestep!" << std::endl;
      }
   }

   void TimestepControlBase::testMaximumTimestep(EPMFloat& rDt)
   {
      // Force maximum timestep
      if(rDt > TimestepConfig::MAXIMUM_TIMESTEP)
      {
         // Set timestep to maximum allowed timestep
         rDt = TimestepConfig::MAXIMUM_TIMESTEP;
      }
   }

   void TimestepControlBase::testCFLCondition(EPMFloat& rDt)
   {
      // Set timestep according to CFL condition
      this->getSimulationCFLCondition();
      if(this->mCFLTimestep > 0.0)
      {
         rDt = std::min(this->mCFLTimestep, rDt);
      }
   }

   void TimestepControlBase::useAdaptiveTimestep(EPMFloat& rDt)
   {
      rDt = this->mController.nextTimestep(this->mError, this->mOldError);
   }

   void TimestepControlBase::useCourantTimestep(EPMFloat& rDt)
   {
      // Include Courant number into timestep value
      rDt = rDt*TimestepConfig::COURANT_NUMBER;
   }

   void TimestepControlBase::useWindowedTimestep(EPMFloat dt)
   {
      // Don't make adaptative scheme to reactive. Allow a range of timestep values
      if(dt < this->rTSParams().dt()*(1.0-TimestepConfig::TIMESTEP_MARGIN) ||  dt > this->rTSParams().dt()*(1.0 + TimestepConfig::TIMESTEP_MARGIN))
      {
         // Store new timestep
         this->rTSParams().updateTimestep(dt);

         std::cout << "Updating Timestep! Now using dt=" << this->rTSParams().dt() << " at time t=" << this->rTSParams().time()<< std::endl;
      }
   }

   void TimestepControlBase::useTimestep(EPMFloat dt)
   {
      // Store new timestep
      this->rTSParams().updateTimestep(dt);

      std::cout << "Updating Timestep! Now using dt=" << this->rTSParams().dt() << " at time t=" << this->rTSParams().time()<< std::endl;
   }

}
