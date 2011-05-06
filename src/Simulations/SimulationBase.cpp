/** \file SimulationBase.cpp
 *  \brief Implementation independent parts of a simulation
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <boost/static_assert.hpp>

// Class include
//
#include "Simulations/SimulationBase.hpp"

// Project includes
//

namespace EPMDynamo {

   SimulationBase::SimulationBase()
      : ComputationBase<IOSystem>(), mSimControl(this->mIOSys.cfg()->aTStep(), this->mEqParams, this->mIOSys.cfg()->aRunI(), this->mIOSys.cfg()->aRun())
   {
   }

   void SimulationBase::initOutput()
   {
      // Write the splitting descriptions
      this->describeSplitting();

      // Initialise all the create writers
      this->mIOSys.initWriters();
   }

   void SimulationBase::writeFiles()
   {
      if(this->mSimControl.tsParams().isNextStep())
      {
         // Save ASCII files
         if(this->mSimControl.tsCounter().triggerASCII())
         {
            this->mIOSys.writeASCII();

            this->mSimControl.checkInterface();
         }

         // Save state file
         if(this->mSimControl.tsCounter().triggerState())
         {
            this->mIOSys.writeHDF5();
         }
      }
   }

   void SimulationBase::preRun()
   {
      // Write loaded state to HDF5 file (this can be important if the truncations were not the same)
      this->mIOSys.writeHDF5();

      // Synchronize simulation
      EPMDYNAMO_SYNCHRONIZE;

      // Execute last minute simulation control initialisation operations
      this->mSimControl.preRun();
   }

   void SimulationBase::postRun()
   {
      // Write the state to HDF5 file before finishing
      this->mIOSys.writeHDF5();

      // Execute early simulation control finalisation operations
      this->mSimControl.postRun();
   }

   void SimulationBase::initFields()
   {
      assert(false); 
   }

   void SimulationBase::addHDF5Output()
   {
      assert(false); 
   }

   void SimulationBase::addASCIIOutput()
   {
      assert(false); 
   }

   void SimulationBase::initEquations()
   {
      assert(false); 
   }

   void SimulationBase::updateEquationsRTP(const int step)
   {
      assert(false); 
   }

   void SimulationBase::updateEquationsRHS()
   {
      assert(false); 
   }

   void SimulationBase::transformEquationsRHS(const int step)
   {
      assert(false); 
   }

   void SimulationBase::addExternalInfluence()
   {
      assert(false); 
   }

   void SimulationBase::timestepEquations()
   {
      assert(false); 
   }

   void SimulationBase::describeSplitting()
   {
      // Create the Physical splitting description file
      PhysicalSplittingFile   physFile(this->mpTrunc);

      // Initialise physical description file
      physFile.init();

      // Write physical description file
      physFile.write();

      // Finalise physical description file
      physFile.finalise();


      // Create the Physical splitting description file
      FDSHSplittingFile   fdshFile(this->mpTrunc);

      // Initialise physical description file
      fdshFile.init();

      // Write physical description file
      fdshFile.write();

      // Finalise physical description file
      fdshFile.finalise();


      // Create the Physical splitting description file
      SpectralSplittingFile   specFile(this->mpTrunc);

      // Initialise physical description file
      specFile.init();

      // Write physical description file
      specFile.write();

      // Finalise physical description file
      specFile.finalise();
   }

}
