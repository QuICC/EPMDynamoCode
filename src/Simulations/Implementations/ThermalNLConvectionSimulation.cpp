/** \file ThermalNLConvectionSimulation.cpp
 *  \brief Implementation of a thermal convection simulation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//

// Project includes
//
#include "Simulations/Implementations/ThermalNLConvectionSimulation.hpp"
#include "BoundaryConditions/L0Harmonic/HeatFluxBC.hpp"

namespace EPMDynamo {

   ThermalNLConvectionSimulation::ThermalNLConvectionSimulation()
      : mCodC(this->mpTrunc, this->mTransform), mVelV(this->mpTrunc, this->mTransform), mTransport(this->mCodC, this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mNavierStokes(this->mVelV, this->mCodC, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void ThermalNLConvectionSimulation::initEquations()
   {
      // Set boundary condition to the transport equation
      // Set constant temperature
      SmartBC  pZeroBC(new ZeroBC(this->mTransform.radBasis()));
      SmartBC  pFluxBC(new HeatFluxBC(-1.0, this->mTransform.radBasis()));
      if(this->mIOSys.cfg()->aBC()(0) == 0)
      {
         this->mTransport.addBC(pZeroBC);
      // Set constant flux
      } else if(this->mIOSys.cfg()->aBC()(0) == 1)
      {
         this->mTransport.addBC(pFluxBC);
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

      // Initialise the transport equation
      this->mTransport.init();

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   void ThermalNLConvectionSimulation::configureTransforms()
   {
      //
      // Setup the SSH transform data manipulator
      //

      // Register transform data packs of transport equation for SSH
      this->registerSSHPacks(this->mTransport.nFSSHPacks(), this->mTransport.nSSHBPacks());

      // Register transform data packs of navier stokes equation for SSH
      this->registerSSHPacks(this->mNavierStokes.nFSSHPacks(), this->mNavierStokes.nSSHBPacks());

      // Configure the SSH manipulator
      this->configureSSHManipulator();

      //
      // Setup the SH transform data manipulator
      //
      
      // Register transform data packs of transport equation for SH
      this->registerSHPacks(this->mTransport.nFSHPacks(), this->mTransport.nSHBPacks());

      // Register transform data packs of navier stokes equation for SH
      this->registerSHPacks(this->mNavierStokes.nFSHPacks(), this->mNavierStokes.nSHBPacks());

      // Configure the SH manipulator
      this->configureSHManipulator();

      //
      // Configure transform nesting
      //

      // Configure the transforms' nesting setup
      this->configureTransformNesting();
   }

   void ThermalNLConvectionSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the transport equation
      this->mTransport.updateRTP(step);

      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   void ThermalNLConvectionSimulation::updateEquationsRHS()
   {
      // Update RHS of the transport equation
      this->mTransport.updateRHS();

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the RTP CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mVelV.oc().rtp());
   }

   void ThermalNLConvectionSimulation::transformEquationsRHS(const int step)
   {
      // Update RHS of the transport equation
      this->mTransport.transformRHS(step);

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   void ThermalNLConvectionSimulation::addExternalInfluence()
   {
      // Add Thermal source term
      this->mTransport.addSourceTerm();
   }

   void ThermalNLConvectionSimulation::timestepEquations()
   {
      // Timestep the transport equation
      this->mTransport.timestep();

      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   void ThermalNLConvectionSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<ThermalNLConvectionTraits> > pInState(new StateFileReader<ThermalNLConvectionTraits>(this->mCodC, this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create a source file reader for the codensity source
      EPMSHARED_PTR<SourceFileReader<ThermalNLConvectionTraits> > pSource(new SourceFileReader<ThermalNLConvectionTraits>(this->mCodC));

      // Read in source state
      this->mIOSys.useSource(pSource);

      // Rescale the source term to correct parameter related value
      this->mCodC.rOcSrc().rescale(this->mEqParams.codSourceScale());

      // Set the kinetic energy scale
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());

   }

   void ThermalNLConvectionSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<ThermalNLConvectionTraits> >  pOutState(new StateFileWriter<ThermalNLConvectionTraits>(this->mCodC, this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   void ThermalNLConvectionSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<EnergyFile<ThermalNLConvectionTraits::CodType> > pCodEnergy(new EnergyFile<ThermalNLConvectionTraits::CodType>(this->mCodC, "cod", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<ThermalNLConvectionTraits::VelType> > pVelEnergy(new EnergyFile<ThermalNLConvectionTraits::VelType>(this->mVelV, "vel", this->mSimControl.tsParams()));
      // Create a energy(m-spectrum) ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<ThermalNLConvectionTraits::VelType> > pVelEnergy_m(new EnergyFile<ThermalNLConvectionTraits::VelType>(this->mVelV, "vel_m", this->mSimControl.tsParams(),2));

      // Add codensity energy to ASCII output
      this->mIOSys.addASCIIWriter(pCodEnergy);
      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);
      this->mIOSys.addASCIIWriter(pVelEnergy_m);

      // Create a energy spectrum ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<SpectrumFile<ThermalNLConvectionTraits::CodType> > pCodSpectrum(new SpectrumFile<ThermalNLConvectionTraits::CodType>(this->mCodC, "cod"));
      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<ThermalNLConvectionTraits::VelType> > pVelSpectrum(new SpectrumFile<ThermalNLConvectionTraits::VelType>(this->mVelV, "vel"));

      // Add codensity energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pCodSpectrum);
      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}
