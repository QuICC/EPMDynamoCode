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
      : mUseASCIIStep((asciiRate > 0)), mUseStateStep((stateRate > 0)), mMaxtstep(maxtstep), mSteps(0), mIsATriggered(false), mIsSTriggered(false), mARate(std::abs(asciiRate)), mSRate(std::abs(stateRate)), mNextAWrite(0), mNextSWrite(0)
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

   void TimestepCounter::increment(const EPMFloat time)
   {
      ++this->mSteps;

      if(this->mUseASCIIStep)
      {
         this->mIsATriggered = (this->mSteps % this->mARate == 0);
      } else
      {
         this->mIsATriggered = (time > this->mNextAWrite);

         if(this->mIsATriggered)
         {
            this->mNextAWrite += this->mARate;
         }
      }

      if(this->mUseStateStep)
      {
         this->mIsSTriggered = (this->mSteps % this->mSRate == 0);
      } else
      {
         this->mIsSTriggered = (time > this->mNextSWrite);

         if(this->mIsSTriggered)
         {
            this->mNextSWrite += this->mSRate;
         }
      }
   }

   void TimestepCounter::printInfo() const
   {
      std::cout << "Simulation stopped after " << this->mSteps << " timesteps." << std::endl;
   }

}
