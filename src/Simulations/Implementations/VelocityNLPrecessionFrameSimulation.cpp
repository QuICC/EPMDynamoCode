/** \file VelocityNLPrecessionFrameSimulation.cpp
 *  \brief Implementation of a velocity diffusion simulation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Simulations/Implementations/VelocityNLPrecessionFrameSimulation.hpp"

// Project includes
//

namespace EPMDynamo {

   VelocityNLPrecessionFrameSimulation::VelocityNLPrecessionFrameSimulation()
      : mVelV(this->mpTrunc, this->mTransform), mNavierStokes(this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void VelocityNLPrecessionFrameSimulation::initEquations()
   {
      // Create Zero boundary condition pointer
      SmartBC  pZeroBC(new ZeroBC(this->mTransform.radBasis()));

      // Set boundary condition to the Navier-Stokes equation
      if(this->mIOSys.cfg()->aBC()(1) == 1)
      {
         SmartBC  pSFBC(new StressFreeTorBC(this->mTransform.radBasis()));
         SmartBC  pDDBC(new DDRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pSFBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDDBC);
      } else if(this->mIOSys.cfg()->aBC()(1) == 2)
      {
         SmartBC  pPrecBC(new PrecessionFrameBC(this->mTransform.radBasis()));
         SmartBC  pDBC(new DRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pPrecBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);
      } else
      {
         SmartBC  pNSBC(new ZeroBC(this->mTransform.radBasis()));
         SmartBC  pDBC(new DRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pNSBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);
      }

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   void VelocityNLPrecessionFrameSimulation::configureTransforms()
   {
      //
      // Setup the SSH transform data manipulator
      //

      // Register transform data packs of navier stokes equation for SSH
      this->registerSSHPacks(this->mNavierStokes.nFSSHPacks(), this->mNavierStokes.nSSHBPacks());

      // Configure the SSH manipulator
      this->configureSSHManipulator();

      //
      // Setup the SH transform data manipulator
      //

      // Register transform data packs of navier stokes equation for SH
      this->registerSHPacks(this->mNavierStokes.nFSHPacks(), this->mNavierStokes.nSHBPacks());

      // Configure the SSH manipulator
      this->configureSHManipulator();

      //
      // Configure transform nesting
      //

      // Configure the transforms' nesting setup
      this->configureTransformNesting();
   }

   void VelocityNLPrecessionFrameSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   void VelocityNLPrecessionFrameSimulation::updateEquationsRHS()
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mVelV.oc().rtp());
   }

   void VelocityNLPrecessionFrameSimulation::transformEquationsRHS(const int step)
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   void VelocityNLPrecessionFrameSimulation::addExternalInfluence()
   {
   }

   void VelocityNLPrecessionFrameSimulation::timestepEquations()
   {
      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   void VelocityNLPrecessionFrameSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<VelocityNLRotDiffusionTraits> > pInState(new StateFileReader<VelocityNLRotDiffusionTraits>(this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Set the energy scale for the magnetic field
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());
   }

   void VelocityNLPrecessionFrameSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<VelocityNLRotDiffusionTraits> >  pOutState(new StateFileWriter<VelocityNLRotDiffusionTraits>(this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   void VelocityNLPrecessionFrameSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<VelocityNLRotDiffusionTraits::VelType> > pVelEnergy(new EnergyFile<VelocityNLRotDiffusionTraits::VelType>(mVelV, "vel", this->mSimControl.tsParams(), 3));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<VelocityNLRotDiffusionTraits::VelType> > pVelSpectrum(new SpectrumFile<VelocityNLRotDiffusionTraits::VelType>(this->mVelV, "vel"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}
