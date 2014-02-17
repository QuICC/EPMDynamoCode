/** \file PrecessionDynamoSimulation.cpp
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
#include "Simulations/Implementations/PrecessionDynamoSimulation.hpp"

// Project includes
//

namespace EPMDynamo {

   PrecessionDynamoSimulation::PrecessionDynamoSimulation()
      : mMagB(this->mpTrunc, this->mTransform), mVelV(this->mpTrunc, this->mTransform), mInduction(this->mMagB, this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mNavierStokes(this->mVelV, this->mMagB, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void PrecessionDynamoSimulation::initEquations()
   {
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

      // Set boundary condition to the induction equation
      if(this->mIOSys.cfg()->aBC()(2) == 0)
      {
         // Set insulator boundary condition
         SmartBC  pInsulatingBC(new InsulatingBC(this->mTransform.radBasis()));

         this->mInduction.addTorBC(pZeroBC);
         this->mInduction.addPolBC(pInsulatingBC);
      } else if(this->mIOSys.cfg()->aBC()(2) == 1)
      {
         // Set  conductor boundary condition
         SmartBC  pConductorTorBC(new ConductorTorBC(this->mTransform.radBasis()));
         SmartBC  pConductorPolBC(new ConductorPolBC(this->mTransform.radBasis()));

         this->mInduction.addTorBC(pConductorTorBC);
         this->mInduction.addPolBC(pConductorPolBC);
      } else
      {
         throw EPMException("PrecessionDynamoSimulation::initEquations", "Did not know what to do with Magnetic BC");
      }

      // Initialise the induction equation
      this->mInduction.init();

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   void PrecessionDynamoSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the induction equation
      this->mInduction.updateRTP(step);

      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   void PrecessionDynamoSimulation::updateEquationsRHS()
   {
      // Update RHS of the induction equation
      this->mInduction.updateRHS();

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mMagB.oc().rtp(), this->mVelV.oc().rtp());
   }

   void PrecessionDynamoSimulation::transformEquationsRHS(const int step)
   {
      // Transform RHS of the induction equation
      this->mInduction.transformRHS(step);

      // Transform RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   void PrecessionDynamoSimulation::addExternalInfluence()
   {
   }

   void PrecessionDynamoSimulation::timestepEquations()
   {
      // Timestep the induction equation
      this->mInduction.timestep();

      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   void PrecessionDynamoSimulation::configureTransforms()
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

   void PrecessionDynamoSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<PrecessionDynamoTraits> > pInState(new StateFileReader<PrecessionDynamoTraits>(this->mMagB, this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Set the energy scale for the magnetic field
      this->mMagB.rOc().rPerturbation().setEnergyScale(this->mEqParams.meFactor());

      // Set the energy scale for the magnetic field
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());
   }

   void PrecessionDynamoSimulation::addHDF5Output()
   {
      // Create a state file writer for regular state output
      EPMSHARED_PTR<StateFileWriter<PrecessionDynamoTraits> >  pOutState(new StateFileWriter<PrecessionDynamoTraits>(this->mMagB, this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      // Add state file writer to hdf5 outputs
      this->mIOSys.addHDF5Writer(pOutState);
   }

   void PrecessionDynamoSimulation::addASCIIOutput()
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
      EPMSHARED_PTR<EnergyFile<PrecessionDynamoTraits::MagType> > pMagEnergy(new EnergyFile<PrecessionDynamoTraits::MagType>(this->mMagB, "mag", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<PrecessionDynamoTraits::VelType> > pVelEnergy(new EnergyFile<PrecessionDynamoTraits::VelType>(this->mVelV, "vel", this->mSimControl.tsParams()));

      // Add magnetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pMagEnergy);
      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<SpectrumFile<PrecessionDynamoTraits::MagType> > pMagSpectrum(new SpectrumFile<PrecessionDynamoTraits::MagType>(this->mMagB, "mag", 1));
      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<PrecessionDynamoTraits::VelType> > pVelSpectrum(new SpectrumFile<PrecessionDynamoTraits::VelType>(this->mVelV, "vel", 1));

      // Add magnetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pMagSpectrum);
      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}
