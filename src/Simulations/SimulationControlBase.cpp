/** \file SimulationControlBase.cpp
 *  \brief Implementation of the simulation control base
 */

// System includes
//

// External includes
//

// Class include
//
#include "Simulations/SimulationControlBase.hpp"

// Project includes
//

namespace EPMDynamo {

   SimulationControlBase::SimulationControlBase(EPMFloat t, EPMFloat dt, int maxtstep, int aRate, int sRate, EPMFloat wall)
      : mKeepRunning(true), mTSParams(t, dt), mTSCounter(maxtstep, aRate, sRate), mRunControl(wall), mUseWall(mTSCounter.unlimited()), mCtrlIface()
   {
   }

   SimulationControlBase::SimulationControlBase(const Array& time, const ArrayI& runI, const Array& run)
      : mKeepRunning(true), mTSParams(time(0), time(1)), mTSCounter(runI(0), runI(1), runI(2)), mRunControl(run(0)), mUseWall(mTSCounter.unlimited()), mCtrlIface()
   {
   }

   void SimulationControlBase::preRun()
   {
      // Set clock start point of run
      this->runControl().start();
   }

   void SimulationControlBase::postRun()
   {
   }

   void SimulationControlBase::printInfo() const
   {
      // Print information from Timestep parameters
      this->tsParams().printInfo();

      // Print information from Timestep counter
      this->tsCounter().printInfo();

      // Print information from Runtime control
      this->runControl().printInfo(this->tsCounter().step());
   }

   void SimulationControlBase::checkInterface()
   {
      // Read input from control interface
      this->mCtrlIface.read();
   }

}
