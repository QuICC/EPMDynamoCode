/** \file LibrationDynamoSimulation.cpp
 *  \brief Implementation of a dynamo simulation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Simulations/Implementations/LibrationDynamoSimulation.hpp"

// Project includes
//

namespace EPMDynamo {

   LibrationDynamoSimulation::LibrationDynamoSimulation()
      : mMagB(this->mpTrunc, this->mTransform), mVelV(this->mpTrunc, this->mTransform), mInduction(this->mMagB, this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mNavierStokes(this->mVelV, this->mMagB, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void LibrationDynamoSimulation::initEquations()
   {
      // Set boundary condition to the Induction equation
      SmartBC  pZeroBC(new ZeroBC(this->mTransform.radBasis()));
      SmartBC  pInsulatingBC(new InsulatingBC(this->mTransform.radBasis()));
      this->mInduction.addTorBC(pZeroBC);
      this->mInduction.addPolBC(pInsulatingBC);

      // Set time dependent longitudinal libration boundary condition
      if(this->mIOSys.cfg()->aBC()(1) == 3)
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
      // Stop if wrong BCs are set for the Navier-Stokes equation
      } else
      {
         throw EPMException("LibrationDynamoSimulation::initEquations", "Simulation requires libration boundary condition!");
      }

      // Initialise the induction equation
      this->mInduction.init();

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   void LibrationDynamoSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the induction equation
      this->mInduction.updateRTP(step);

      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   void LibrationDynamoSimulation::updateEquationsRHS()
   {
      // Update RHS of the induction equation
      this->mInduction.updateRHS();

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mMagB.oc().rtp(), this->mVelV.oc().rtp());
   }

   void LibrationDynamoSimulation::transformEquationsRHS(const int step)
   {
      // Transform RHS of the induction equation
      this->mInduction.transformRHS(step);

      // Transform RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   void LibrationDynamoSimulation::addExternalInfluence()
   {
   }

   void LibrationDynamoSimulation::timestepEquations()
   {
      // Timestep the induction equation
      this->mInduction.timestep();

      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   void LibrationDynamoSimulation::configureTransforms()
   {
      //
      // Setup the SSH transform data manipulator
      //
      
      // Register transform data packs of induction equation for SSH
      this->registerSSHPacks(this->mInduction.nFSSHPacks(), this->mInduction.nSSHBPacks());

      // Register transform data packs of navier stokes equation for SSH
      this->registerSSHPacks(this->mNavierStokes.nFSSHPacks(), this->mNavierStokes.nSSHBPacks());

      // Configure the SSH manipulator
      this->configureSSHManipulator();

      //
      // Setup the SH transform data manipulator
      //
      
      // Register transform data packs of induction equation for SH
      this->registerSHPacks(this->mInduction.nFSHPacks(), this->mInduction.nSHBPacks());

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

   void LibrationDynamoSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<LibrationDynamoTraits> > pInState(new StateFileReader<LibrationDynamoTraits>(this->mMagB, this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Set the energy scale for the magnetic field
      this->mMagB.rOc().rPerturbation().setEnergyScale(this->mEqParams.meFactor());

      // Set the energy scale for the magnetic field
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());
   }

   void LibrationDynamoSimulation::addHDF5Output()
   {
      // Create a state file writer for regular state output
      EPMSHARED_PTR<StateFileWriter<LibrationDynamoTraits> >  pOutState(new StateFileWriter<LibrationDynamoTraits>(this->mMagB, this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      // Add state file writer to hdf5 outputs
      this->mIOSys.addHDF5Writer(pOutState);
   }

   void LibrationDynamoSimulation::addASCIIOutput()
   {
      // Create a CFL timestep ASCII logging file
      EPMSHARED_PTR<CFLTimestepFile> pCFLFile(new CFLTimestepFile("cfl_timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pCFLFile);

      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<EnergyFile<LibrationDynamoTraits::MagType> > pMagEnergy(new EnergyFile<LibrationDynamoTraits::MagType>(this->mMagB, "mag", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<LibrationDynamoTraits::VelType> > pVelEnergy(new EnergyFile<LibrationDynamoTraits::VelType>(this->mVelV, "vel", this->mSimControl.tsParams()));

      // Add magnetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pMagEnergy);
      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<SpectrumFile<LibrationDynamoTraits::MagType> > pMagSpectrum(new SpectrumFile<LibrationDynamoTraits::MagType>(this->mMagB, "mag"));
      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<LibrationDynamoTraits::VelType> > pVelSpectrum(new SpectrumFile<LibrationDynamoTraits::VelType>(this->mVelV, "vel"));

      // Add magnetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pMagSpectrum);
      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}
