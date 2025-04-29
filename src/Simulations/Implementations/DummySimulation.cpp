/** \file DummySimulation.cpp
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
#include "Simulations/Implementations/DummySimulation.hpp"

// Project includes
//
#include "BoundaryConditions/L0Harmonic/HeatFluxBC.hpp"

namespace EPMDynamo {

   DummySimulation::DummySimulation()
      : mCodC(this->mpTrunc, this->mTransform), mMagB(this->mpTrunc, this->mTransform), mVelV(this->mpTrunc, this->mTransform), mInduction(this->mMagB, this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mTransport(this->mCodC, this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mNavierStokes(this->mVelV, this->mMagB, this->mCodC, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void DummySimulation::initEquations()
   {
      SmartBC  pZeroBC(new ZeroBC(this->mTransform.radBasis()));

      // Set boundary condition to the transport equation
      if(this->mIOSys.cfg()->aBC()(0) == 0)
      {
         // Set constant temperature
         this->mTransport.addBC(pZeroBC);
      // Set constant flux
      } else if(this->mIOSys.cfg()->aBC()(0) == 1)
      {
         SmartBC  pFluxBC(new HeatFluxBC(-1.0, this->mTransform.radBasis()));

         this->mTransport.addBC(pFluxBC);
      } else
      {
         throw EPMException("DummySimulation::initEquations", "Did not know what to do with Codensity BC");
      }

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
      } else if(this->mIOSys.cfg()->aBC()(1) == 0)
      {
         SmartBC  pNSBC(new ZeroBC(this->mTransform.radBasis()));
         SmartBC  pDBC(new DRadialBC(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pNSBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);
      } else
      {
         throw EPMException("DummySimulation::initEquations", "Did not know what to do with Velocity BC");
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
         throw EPMException("DummySimulation::initEquations", "Did not know what to do with Magnetic BC");
      }

      // Initialise the induction equation
      this->mInduction.init();

      // Initialise the transport equation
      this->mTransport.init();

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   void DummySimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the induction equation
      this->mInduction.updateRTP(step);

      // Update RTP values of the transport equation
      this->mTransport.updateRTP(step);

      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   void DummySimulation::updateEquationsRHS()
   {
      // Update RHS of the induction equation
      this->mInduction.updateRHS();

      // Update RHS of the transport equation
      this->mTransport.updateRHS();

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mMagB.oc().rtp(), this->mVelV.oc().rtp());
   }

   void DummySimulation::transformEquationsRHS(const int step)
   {
      // Transform RHS of the induction equation
      this->mInduction.transformRHS(step);

      // Transform RHS of the transport equation
      this->mTransport.transformRHS(step);

      // Transform RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   void DummySimulation::addExternalInfluence()
   {
      // Add Thermal source term
      this->mTransport.addSourceTerm();
   }

   void DummySimulation::timestepEquations()
   {
      // Timestep the induction equation
      this->mInduction.timestep();

      // Timestep the transport equation
      this->mTransport.timestep();

      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   void DummySimulation::configureTransforms()
   {
      //
      // Setup the SSH transform data manipulator
      //
      
      // Register transform data packs of induction equation for SSH
      this->registerSSHPacks(this->mInduction.nFSSHPacks(), this->mInduction.nSSHBPacks());

      // Register transform data packs of transport equation for SSH
      this->registerSSHPacks(this->mTransport.nFSSHPacks(), this->mTransport.nSSHBPacks());

      // Register transform data packs of navier stokes equation for SSH
      this->registerSSHPacks(this->mNavierStokes.nFSSHPacks(), this->mNavierStokes.nSSHBPacks());

      // Configure the SSH manipulator
      this->configureSSHManipulator();

      //
      // Setup the SH transform data manipulator
      //
      
      // Register transform data packs of induction equation for SH
      this->registerSHPacks(this->mInduction.nFSHPacks(), this->mInduction.nSHBPacks());

      // Register transform data packs of transport equation for SH
      this->registerSHPacks(this->mTransport.nFSHPacks(), this->mTransport.nSHBPacks());

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

   void DummySimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<DynamoTraits> > pInState(new StateFileReader<DynamoTraits>(this->mCodC, this->mMagB, this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create a source file reader for the codensity source
      EPMSHARED_PTR<SourceFileReader<DynamoTraits> > pSource(new SourceFileReader<DynamoTraits>(this->mCodC));

      // Read in source state
      this->mIOSys.useSource(pSource);

      // Rescale source term to correct values for given parameters
      this->mCodC.rOcSrc().rescale(this->mEqParams.codSourceScale());

      // Set the energy scale for the magnetic field
      this->mMagB.rOc().rPerturbation().setEnergyScale(this->mEqParams.meFactor());

      // Set the energy scale for the magnetic field
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());
   }

   void DummySimulation::addHDF5Output()
   {
      // Create a state file writer for regular state output
      EPMSHARED_PTR<StateFileWriter<DynamoTraits> >  pOutState(new StateFileWriter<DynamoTraits>(this->mCodC, this->mMagB, this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      // Add state file writer to hdf5 outputs
      this->mIOSys.addHDF5Writer(pOutState);
   }

   void DummySimulation::addASCIIOutput()
   {
      // Create a CFL timestep ASCII logging file
      EPMSHARED_PTR<CFLTimestepFile> pCFLFile(new CFLTimestepFile("cfl_timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pCFLFile);

      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<EnergyFile<DynamoTraits::CodType> > pCodEnergy(new EnergyFile<DynamoTraits::CodType>(this->mCodC, "cod", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<EnergyFile<DynamoTraits::MagType> > pMagEnergy(new EnergyFile<DynamoTraits::MagType>(this->mMagB, "mag", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<DynamoTraits::VelType> > pVelEnergy(new EnergyFile<DynamoTraits::VelType>(this->mVelV, "vel", this->mSimControl.tsParams()));

      // Add codensity energy to ASCII output
      this->mIOSys.addASCIIWriter(pCodEnergy);
      // Add magnetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pMagEnergy);
      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<SpectrumFile<DynamoTraits::CodType> > pCodSpectrum(new SpectrumFile<DynamoTraits::CodType>(this->mCodC, "cod"));
      // Create a energy spectrum ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<SpectrumFile<DynamoTraits::MagType> > pMagSpectrum(new SpectrumFile<DynamoTraits::MagType>(this->mMagB, "mag"));
      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<DynamoTraits::VelType> > pVelSpectrum(new SpectrumFile<DynamoTraits::VelType>(this->mVelV, "vel"));

      // Add codensity energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pCodSpectrum);
      // Add magnetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pMagSpectrum);
      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}
