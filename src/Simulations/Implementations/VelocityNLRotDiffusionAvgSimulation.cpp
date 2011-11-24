/** \file VelocityNLRotDiffusionAvgSimulation.cpp
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
#include "Simulations/Implementations/VelocityNLRotDiffusionAvgSimulation.hpp"

// Project includes
//
#include "IO/ASCII/ExtractRadiusFile.hpp"

namespace EPMDynamo {

   VelocityNLRotDiffusionAvgSimulation::VelocityNLRotDiffusionAvgSimulation()
      : mVelV(this->mpTrunc, this->mTransform), mNavierStokes(this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mTimeAverager(this->mVelV, this->mTransform, this->mSimControl.tsParams())
   {
   }

   void VelocityNLRotDiffusionAvgSimulation::initEquations()
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

      // Set precession frame boundary condition
      } else if(this->mIOSys.cfg()->aBC()(1) == 2)
      {
         SmartBC  pPrecBC(new PrecessionFrameBC(this->mTransform.radBasis()));
         SmartBC  pDBC(new DRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pPrecBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);

      // Set time dependent longitudinal libration boundary condition
      } else if(this->mIOSys.cfg()->aBC()(1) == 3)
      {
         SmartBC  pLibBC(new LongitudinalLibrationBC(this->mTransform.radBasis(), this->mSimControl.tsParams()));
         SmartBC  pDBC(new DRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pLibBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);

      // Set time dependent latitudinal libration boundary condition
      } else if(this->mIOSys.cfg()->aBC()(1) == 4)
      {
         SmartBC  pLibBC(new LatitudinalLibrationBC(this->mTransform.radBasis(), this->mSimControl.tsParams()));
         SmartBC  pDBC(new DRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pLibBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);
      } else
      {
         throw EPMException("VelocityNLRotDiffusionSimuation::initEquations", "Did not know what to do with Velocity BC");
      }

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();

      // Initialise the TimeAverager
      this->mTimeAverager.init();
   }

   void VelocityNLRotDiffusionAvgSimulation::configureTransforms()
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

   void VelocityNLRotDiffusionAvgSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   void VelocityNLRotDiffusionAvgSimulation::updateEquationsRHS()
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mVelV.oc().rtp());
   }

   void VelocityNLRotDiffusionAvgSimulation::transformEquationsRHS(const int step)
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   void VelocityNLRotDiffusionAvgSimulation::addExternalInfluence()
   {
   }

   void VelocityNLRotDiffusionAvgSimulation::timestepEquations()
   {
      // Update the time average computation
      this->mTimeAverager.timestep();

      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   void VelocityNLRotDiffusionAvgSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<VelocityNLRotDiffusionTraits> > pInState(new StateFileReader<VelocityNLRotDiffusionTraits>(this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Set the energy scale for the magnetic field
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());
   }

   void VelocityNLRotDiffusionAvgSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<VelocityNLRotDiffusionTraits> >  pOutState(new StateFileWriter<VelocityNLRotDiffusionTraits>(this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);

      EPMSHARED_PTR<StateFileWriter<VelocityNLRotDiffusionTraits> >  pAvgState(new StateFileWriter<VelocityNLRotDiffusionTraits>(*this->mTimeAverager.pAvgVelV(), this->mEqParams, this->mSimControl.tsParams(), "Avg"));

      this->mIOSys.addHDF5Writer(pAvgState);
   }

   void VelocityNLRotDiffusionAvgSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<VelocityNLRotDiffusionTraits::VelType> > pVelEnergy(new EnergyFile<VelocityNLRotDiffusionTraits::VelType>(mVelV, "vel", this->mSimControl.tsParams()));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy ASCII diagnostic file for the averaged velocity field
      EPMSHARED_PTR<EnergyFile<VelocityNLRotDiffusionTraits::VelType> > pAvgVelEnergy(new EnergyFile<VelocityNLRotDiffusionTraits::VelType>(*this->mTimeAverager.pAvgVelV(), "avgVel", this->mSimControl.tsParams()));

      // Add averaged kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pAvgVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<VelocityNLRotDiffusionTraits::VelType> > pVelSpectrum(new SpectrumFile<VelocityNLRotDiffusionTraits::VelType>(this->mVelV, "vel", 1));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);

      // Create a averaged energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<VelocityNLRotDiffusionTraits::VelType> > pAvgVelSpectrum(new SpectrumFile<VelocityNLRotDiffusionTraits::VelType>(*this->mTimeAverager.pAvgVelV(), "avgVel"));

      // Add averaged kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pAvgVelSpectrum);

      // Add libration output file if required
      if(this->mIOSys.cfg()->aBC()(1) == 4)
      {
         // Create a libration ASCII diagnostic file for the velocity field
         EPMSHARED_PTR<LibrationFile<VelocityNLRotDiffusionTraits::VelType> > pLibration(new LibrationFile<VelocityNLRotDiffusionTraits::VelType>(this->mVelV, "vel", this->mSimControl.tsParams(), this->mTransform.radBasis()));

         // Add kinetic energy spectrum to ASCII output
         this->mIOSys.addASCIIWriter(pLibration);
      }
   }

}
