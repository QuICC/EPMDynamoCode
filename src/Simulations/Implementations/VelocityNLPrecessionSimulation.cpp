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
#include "IO/ASCII/ExtractRadiusFile.hpp"

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
      if(this->mIOSys.cfg()->aBC()(1) != 0)
      {
         throw EPMException("PrecessionDynamoSimulation::initEquations", "Simulation requires no-slip boundary condition!");
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
      EPMSHARED_PTR<StateFileReader<VelocityNLRotDiffusionTraits> > pInState(new StateFileReader<VelocityNLRotDiffusionTraits>(this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Set the energy scale for the magnetic field
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());
   }

   void VelocityNLPrecessionSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<VelocityNLRotDiffusionTraits> >  pOutState(new StateFileWriter<VelocityNLRotDiffusionTraits>(this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   void VelocityNLPrecessionSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<VelocityNLRotDiffusionTraits::VelType> > pVelEnergy(new EnergyFile<VelocityNLRotDiffusionTraits::VelType>(mVelV, "vel", this->mSimControl.tsParams(), 3));

      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<VelocityNLRotDiffusionTraits::VelType> > pVelEnergy_M(new EnergyFile<VelocityNLRotDiffusionTraits::VelType>(mVelV, "vel_m", this->mSimControl.tsParams(), 2));
      
      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);
      this->mIOSys.addASCIIWriter(pVelEnergy_M);

      // Create a solib body rotation ASCII diagnostic file for the velocity field
         EPMSHARED_PTR<LibrationFile<VelocityNLRotDiffusionTraits::VelType> > pLibration(new LibrationFile<VelocityNLRotDiffusionTraits::VelType>(this->mVelV, "vel", this->mSimControl.tsParams(), this->mTransform.radBasis()));

      // Add solid body rotation to ASCII output
         this->mIOSys.addASCIIWriter(pLibration);



      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<VelocityNLRotDiffusionTraits::VelType> > pVelSpectrum(new SpectrumFile<VelocityNLRotDiffusionTraits::VelType>(this->mVelV, "vel",1));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);

      // Add velocity probe
      Array ratios(3);
      ratios(0) = 0.4;
      ratios(1) = 1.0/2.0;
      ratios(2) = 0.0;
      // Creat velocity probe file
      EPMSHARED_PTR<FieldProbeFile<VelocityNLRotDiffusionTraits::VelType> > pVelProbe1(new FieldProbeFile<VelocityNLRotDiffusionTraits::VelType>(this->mVelV, "vel1", this->mSimControl.tsParams(), ratios));
      // Add velocity probe to ASCII output
      this->mIOSys.addASCIIWriter(pVelProbe1);

      ratios(0) = 0.6;
      ratios(1) = 1.0/2.0;
      ratios(2) = 0.0;
      // Creat velocity probe file
      EPMSHARED_PTR<FieldProbeFile<VelocityNLRotDiffusionTraits::VelType> > pVelProbe2(new FieldProbeFile<VelocityNLRotDiffusionTraits::VelType>(this->mVelV, "vel2", this->mSimControl.tsParams(), ratios));
      // Add velocity probe to ASCII output
      this->mIOSys.addASCIIWriter(pVelProbe2);

      
      ratios(0) = 0.9;
      ratios(1) = 1.0/2.0;
      ratios(2) = 0.0;
      // Creat velocity probe file
      EPMSHARED_PTR<FieldProbeFile<VelocityNLRotDiffusionTraits::VelType> > pVelProbe3(new FieldProbeFile<VelocityNLRotDiffusionTraits::VelType>(this->mVelV, "vel3", this->mSimControl.tsParams(), ratios));
      // Add velocity probe to ASCII output
      this->mIOSys.addASCIIWriter(pVelProbe3);
   }
}
