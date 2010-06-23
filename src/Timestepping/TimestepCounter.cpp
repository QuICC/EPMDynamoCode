/** \file TimestepCounter.cpp
 *  \brief Source of the implementation of the timestep data counter
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/TimestepCounter.hpp"

// Project includes
//

namespace EPMDynamo {

   TimestepCounter::TimestepCounter(int maxtstep, int asciiRate, int stateRate)
      : mMaxtstep(maxtstep), mSteps(0), mARate(asciiRate), mSRate(stateRate)
   {
   }

   bool TimestepCounter::keepRunning() const
   {
      if(this->mSteps >= this->mMaxtstep)
      {
         std::cout << "Maxtstep reached!" <<  std::endl;
         return false;
      } else
      {
         return true;
      }
   }

   bool TimestepCounter::isUnlimited() const
   {
      if(this->mMaxtstep < 1)
      {
         return true;
      } else
      {
         return false;
      }
   }

   bool TimestepCounter::triggerASCII() const
   {
      if(this->mSteps % this->mARate == 0)
      {
         return true;
      } else
      {
         return false;
      }
   }

   bool TimestepCounter::triggerState() const
   {
      if(this->mSteps % this->mSRate == 0)
      {
         return true;
      } else
      {
         return false;
      }
   }

   void TimestepCounter::printInfo() const
   {
      std::cout << "Simulation stopped after " << this->mSteps << " timesteps." << std::endl;
   }

}
