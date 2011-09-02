/** \file VelocityRotDiffusionSimulation.cpp
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
#include "Simulations/Implementations/VelocityRotDiffusionSimulation.hpp"

// Project includes
//

namespace EPMDynamo {

   VelocityRotDiffusionSimulation::VelocityRotDiffusionSimulation()
      : mVelV(this->mpTrunc, this->mTransform), mNavierStokes(this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void VelocityRotDiffusionSimulation::initEquations()
   {
      // Create Zero boundary condition pointer
      SmartBC  pZeroBC(new ZeroBC(this->mTransform.radBasis()));

      // Set boundary condition to the Navier-Stokes equation
      //

      // Set no-slip boundary condition
      if(this->mIOSys.cfg()->aBC()(1) == 0)
      {
         SmartBC  pNSBC(new ZeroBC(this->mTransform.radBasis()));
         SmartBC  pDBC(new DRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pNSBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);

      // Set stress-free boundary condition
      } else if(this->mIOSys.cfg()->aBC()(1) == 1)
      {
         SmartBC  pSFBC(new StressFreeTorBC(this->mTransform.radBasis()));
         SmartBC  pDDBC(new DDRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pSFBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDDBC);

      // Set time dependent longitudinal libration boundary condition
      } else if(this->mIOSys.cfg()->aBC()(1) == 3)
      {
         SmartBC  pPrecBC(new LongitudinalLibrationBC(this->mTransform.radBasis(), this->mSimControl.tsParams()));
         SmartBC  pDBC(new DRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pPrecBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);

      // Set time dependent latitudinal libration boundary condition
      } else if(this->mIOSys.cfg()->aBC()(1) == 4)
      {
         SmartBC  pPrecBC(new LatitudinalLibrationBC(this->mTransform.radBasis(), this->mSimControl.tsParams()));
         SmartBC  pDBC(new DRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pPrecBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);
      } else
      {
         throw EPMException("VelocityRotDiffusionSimuation::initEquations", "Did not know what to do with Velocity BC");
      }

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   void VelocityRotDiffusionSimulation::configureTransforms()
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

   void VelocityRotDiffusionSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   void VelocityRotDiffusionSimulation::updateEquationsRHS()
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mVelV.oc().rtp());
   }

   void VelocityRotDiffusionSimulation::transformEquationsRHS(const int step)
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   void VelocityRotDiffusionSimulation::addExternalInfluence()
   {
   }

   void VelocityRotDiffusionSimulation::timestepEquations()
   {
      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   void VelocityRotDiffusionSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<VelocityRotDiffusionTraits> > pInState(new StateFileReader<VelocityRotDiffusionTraits>(this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Set the energy scale for the magnetic field
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());
   }

   void VelocityRotDiffusionSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<VelocityRotDiffusionTraits> >  pOutState(new StateFileWriter<VelocityRotDiffusionTraits>(this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   void VelocityRotDiffusionSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<VelocityRotDiffusionTraits::VelType> > pVelEnergy(new EnergyFile<VelocityRotDiffusionTraits::VelType>(mVelV, "vel", this->mSimControl.tsParams()));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<VelocityRotDiffusionTraits::VelType> > pVelSpectrum(new SpectrumFile<VelocityRotDiffusionTraits::VelType>(this->mVelV, "vel"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);

      // Add libration output file if required
      if(this->mIOSys.cfg()->aBC()(1) == 4)
      {
         // Create a libration ASCII diagnostic file for the velocity field
         EPMSHARED_PTR<LibrationFile<VelocityRotDiffusionTraits::VelType> > pLibration(new LibrationFile<VelocityRotDiffusionTraits::VelType>(this->mVelV, "vel", this->mSimControl.tsParams(), this->mTransform.radBasis()));

         // Add kinetic energy spectrum to ASCII output
         this->mIOSys.addASCIIWriter(pLibration);
      }
   }

}
