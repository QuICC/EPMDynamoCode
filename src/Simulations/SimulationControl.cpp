/** \file SimulationControl.cpp
 *  \brief Implementation of a general simulation control structure
 */

// System includes
//

// External includes
//

// Class include
//
#include "Simulations/SimulationControl.hpp"

// Project includes
//

namespace EPMDynamo {

   SimulationControl::SimulationControl(EPMFloat t, EPMFloat dt, const EquationParameters& eqParams, int maxtstep, int aRate, int sRate, EPMFloat wall, EPMFloat maxtime)
      : SimulationControlBase(t, dt, maxtstep, aRate, sRate, wall, maxtime), mTSControl(mTSParams, eqParams)
   {
   }

   SimulationControl::SimulationControl(const Array& time, const EquationParameters& eqParams, const ArrayI& runI, const Array& run)
      : SimulationControlBase(time, runI, run), mTSControl(mTSParams, eqParams)
   {
   }

   void SimulationControl::update()
   {
      // Check convergence of time integration
      this->tsControl().checkConvergence(this->tsCounter().steps());

      if(this->tsParams().isNextStep())
      {
         // Update simulation time
         this->tsParams().updateTime();

         // increment the number of performed timesteps
         this->tsCounter().increment();

         // Update runtime
         this->runControl().update(this->tsParams().time());
      }

      // Keep running simulation?
      if(this->mUseWall)
      {
         this->mKeepRunning = this->tsControl().keepRunning() && this->runControl().keepRunning() && this->ctrlIface().keepRunning();
      } else
      {
         this->mKeepRunning = this->tsControl().keepRunning() && this->tsCounter().keepRunning() && this->ctrlIface().keepRunning();
      }
   }
}
