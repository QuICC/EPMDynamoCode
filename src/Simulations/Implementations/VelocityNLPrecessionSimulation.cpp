/** \file VelocityNLPrecessionSimulation.cpp
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
#include "Simulations/Implementations/VelocityNLPrecessionSimulation.hpp"

// Project includes
//

namespace EPMDynamo {

   VelocityNLPrecessionSimulation::VelocityNLPrecessionSimulation()
      : mVelV(this->mpTrunc, this->mTransform), mNavierStokes(this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void VelocityNLPrecessionSimulation::initEquations()
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
         SmartBC  pPrecBC(new PrecessionBC(this->mTransform.radBasis()));
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

   void VelocityNLPrecessionSimulation::configureTransforms()
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

   void VelocityNLPrecessionSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   void VelocityNLPrecessionSimulation::updateEquationsRHS()
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mVelV.oc().rtp());
   }

   void VelocityNLPrecessionSimulation::transformEquationsRHS(const int step)
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   void VelocityNLPrecessionSimulation::addExternalInfluence()
   {
   }

   void VelocityNLPrecessionSimulation::timestepEquations()
   {
      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   void VelocityNLPrecessionSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<VelocityNLPrecessionTraits> > pInState(new StateFileReader<VelocityNLPrecessionTraits>(this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());
   }

   void VelocityNLPrecessionSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<VelocityNLPrecessionTraits> >  pOutState(new StateFileWriter<VelocityNLPrecessionTraits>(this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   void VelocityNLPrecessionSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<VelocityNLPrecessionTraits::VelType> > pVelEnergy(new EnergyFile<VelocityNLPrecessionTraits::VelType>(mVelV, "vel", this->mSimControl.tsParams()));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<VelocityNLPrecessionTraits::VelType> > pVelSpectrum(new SpectrumFile<VelocityNLPrecessionTraits::VelType>(this->mVelV, "vel"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}
