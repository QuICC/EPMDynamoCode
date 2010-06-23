/** \file TimestepParameters.cpp
 *  \brief Source of the implementation of the timestep parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/TimestepParameters.hpp"

// Project includes
//
#include "Timestepping/TimestepConfig.hpp"

namespace EPMDynamo {

   TimestepParameters::TimestepParameters(EPMFloat time, EPMFloat timestep)
      : mIsNextStep(true), mHasNewDt(true), mIsRejected(false), mTime(time), mTimestep(timestep), mOldTimestep(timestep), mError(0.0)
   {
      // Check initial values
      this->checkValues();
   }
 
   void TimestepParameters::updateError(EPMFloat err)
   {
      this->mError = err;
   }
 
   void TimestepParameters::resetError()
   {
      this->mError = 0.0;
   }

   void TimestepParameters::checkValues()
   {
      // Check for initiatisation values
      if(this->mTimestep <= 0.0)
      {
         this->mTimestep = TimestepConfig::MIN_TIMESTEP;
         this->mOldTimestep = TimestepConfig::MIN_TIMESTEP;
      }
      if(this->mTime < 0.0)
      {
         this->mTime = 0.0;
      }
   }

   void TimestepParameters::updateTimestep(EPMFloat dt)
   {
      // Store previous timestep
      this->mOldTimestep = this->mTimestep;

      // Set new timestep
      this->mTimestep = dt;

      // Trigger timestep matrices update
      this->mHasNewDt = true;
   }

   void TimestepParameters::printInfo() const
   {
      std::cout << "Simulation time reached is: " << this->mTime << std::endl;
   }

}
