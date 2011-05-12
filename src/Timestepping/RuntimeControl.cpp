/** \file RuntimeControl.cpp
 *  \brief Source of the implementation of a simple run time control
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "Timestepping/RuntimeControl.hpp"

// Project includes
//

namespace EPMDynamo {

   RuntimeControl::RuntimeControl(EPMFloat wall, EPMFloat maxtime)
      : mWallTime(wall), mMaxtime(maxtime), mRuntime(0.0), mTime(0.0)
   {
   }

   bool RuntimeControl::keepRunning() const
   {
      if(this->mWallTime > 0 && this->mRuntime >= this->mWallTime)
      {
         std::cout << "Wall time reached!" << std::endl;

         return false;
      } else if(this->mMaxtime > 0 && this->mTime >= this->mMaxtime)
      {
         std::cout << "Max integration time reached!" << std::endl;

         return false;
      } else
      {
         return true;
      }
   }

   void RuntimeControl::update(EPMFloat time)
   {
      this->mRuntime += this->resetTimer()/3600.0;

      this->mTime = time;
   }

   void RuntimeControl::printInfo(const int steps) const
   {
      std::cout << "Total runtime was: "<< this->mRuntime << " hours (" << this->mRuntime*3600 << " seconds)" << std::endl;
      std::cout << "Average timestep time was: " << this->mRuntime*3600. / static_cast<EPMFloat>(steps) << " seconds/step" << std::endl;
   }

}
