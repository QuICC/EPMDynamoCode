/** \file MagneticDiffusionSimulation.cpp
 *  \brief Implementation of a magnetic diffusion simulation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Simulations/Implementations/MagneticDiffusionSimulation.hpp"

// Project includes
//

namespace EPMDynamo {

   MagneticDiffusionSimulation::MagneticDiffusionSimulation()
      : mMagB(this->mpTrunc, this->mTransform), mInduction(mMagB, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void MagneticDiffusionSimulation::initEquations()
   {
      // Set boundary condition to the Induction equation
      SmartBC  pZeroBC(new ZeroBC(this->mTransform.radBasis()));
      SmartBC  pInsulatingBC(new InsulatingBC(this->mTransform.radBasis()));
      this->mInduction.addTorBC(pZeroBC);
      this->mInduction.addPolBC(pInsulatingBC);

      // Initialise the Induction equation
      this->mInduction.init();
   }

   void MagneticDiffusionSimulation::configureTransforms()
   {
      //
      // Setup the SSH transform data manipulator
      //

      // Register transform data packs of induction equation for SSH
      this->registerSSHPacks(this->mInduction.nFSSHPacks(), this->mInduction.nSSHBPacks());

      // Configure the SSH manipulator
      this->configureSSHManipulator();

      //
      // Setup the SH transform data manipulator
      //

      // Register transform data packs of induction equation for SH
      this->registerSHPacks(this->mInduction.nFSHPacks(), this->mInduction.nSHBPacks());

      // Configure the SSH manipulator
      this->configureSHManipulator();

      //
      // Configure transform nesting
      //

      // Configure the transforms' nesting setup
      this->configureTransformNesting();
   }

   void MagneticDiffusionSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the Induction equation
      this->mInduction.updateRTP(step);
   }

   void MagneticDiffusionSimulation::updateEquationsRHS()
   {
      // Update RHS of the induction equation
      this->mInduction.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mMagB.oc().rtp());
   }

   void MagneticDiffusionSimulation::transformEquationsRHS(const int step)
   {
      // Update RHS of the induction equation
      this->mInduction.transformRHS(step);
   }

   void MagneticDiffusionSimulation::addExternalInfluence()
   {
   }

   void MagneticDiffusionSimulation::timestepEquations()
   {
      // Timestep the induction equation
      this->mInduction.timestep();
   }

   void MagneticDiffusionSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<MagneticDiffusionTraits> > pInState(new StateFileReader<MagneticDiffusionTraits>(this->mMagB,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());
   }

   void MagneticDiffusionSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<MagneticDiffusionTraits> >  pOutState(new StateFileWriter<MagneticDiffusionTraits>(this->mMagB, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   void MagneticDiffusionSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<EnergyFile<MagneticDiffusionTraits::MagType> > pMagEnergy(new EnergyFile<MagneticDiffusionTraits::MagType>(mMagB, "mag", this->mSimControl.tsParams()));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pMagEnergy);

      // Create a energy spectrum ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<SpectrumFile<MagneticDiffusionTraits::MagType> > pMagSpectrum(new SpectrumFile<MagneticDiffusionTraits::MagType>(this->mMagB, "mag"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pMagSpectrum);
   }

}
