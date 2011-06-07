/** \file ChandrasekharOnsetSimulation.cpp
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
#include "Simulations/Implementations/ChandrasekharOnsetSimulation.hpp"

// Project includes
//

namespace EPMDynamo {

   ChandrasekharOnsetSimulation::ChandrasekharOnsetSimulation()
      : mCodC(this->mpTrunc, this->mTransform), mVelV(this->mpTrunc, this->mTransform), mTransport(this->mCodC, this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mNavierStokes(this->mVelV, this->mCodC, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void ChandrasekharOnsetSimulation::initEquations()
   {
      // Set boundary condition to the transport equation
      SmartBC  pZeroBC(new ZeroBC(this->mTransform.radBasis()));
      this->mTransport.addBC(pZeroBC);

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

   void ChandrasekharOnsetSimulation::configureTransforms()
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

   void ChandrasekharOnsetSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the transport equation
      this->mTransport.updateRTP(step);

      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   void ChandrasekharOnsetSimulation::updateEquationsRHS()
   {
      // Update RHS of the transport equation
      this->mTransport.updateRHS();

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the RTP CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mVelV.oc().rtp());
   }

   void ChandrasekharOnsetSimulation::transformEquationsRHS(const int step)
   {
      // Update RHS of the transport equation
      this->mTransport.transformRHS(step);

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   void ChandrasekharOnsetSimulation::addExternalInfluence()
   {
      // Add Thermal source term
      this->mTransport.addSourceTerm();
   }

   void ChandrasekharOnsetSimulation::timestepEquations()
   {
      // Timestep the transport equation
      this->mTransport.timestep();

      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();

      // Restrict the fields to a specific (set of ) harmonic degree
      this->restrictSolution();
   }

   void ChandrasekharOnsetSimulation::restrictSolution()
   {
      // Lowest harmonic degree to keep in computation
      int l0 = 15;

      // Get the core store harmonic degrees
      ArrayI ls = this->mVelV.oc().trunc()->local()->spec()->lArray();

      // Loop over all the harmonic degrees
      for(int i = 0; i < ls.size(); ++i)
      {
         if(ls(i) < l0)
         {
            // Set Codensity values to zero
            this->mCodC.rOc().rPerturbation().rLShell(i).setConstant(EPMComplex(0.0,0.0));

            // Set Velocity values to zero
            this->mVelV.rOc().rPerturbation().rTor().rLShell(i).setConstant(EPMComplex(0.0,0.0));
            this->mVelV.rOc().rPerturbation().rPol().rLShell(i).setConstant(EPMComplex(0.0,0.0));
         }
      }
   }

   void ChandrasekharOnsetSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<LinearThermalConvectionTraits> > pInState(new StateFileReader<LinearThermalConvectionTraits>(this->mCodC, this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create a source file reader for the codensity source
      EPMSHARED_PTR<SourceFileReader<LinearThermalConvectionTraits> > pSource(new SourceFileReader<LinearThermalConvectionTraits>(this->mCodC));

      // Read in source state
      this->mIOSys.useSource(pSource);

      // Rescale the source term to correct parameter related value
      this->mCodC.rOcSrc().rescale(this->mEqParams.codSourceScale());

      // Set the kinetic energy scale
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());

   }

   void ChandrasekharOnsetSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<LinearThermalConvectionTraits> >  pOutState(new StateFileWriter<LinearThermalConvectionTraits>(this->mCodC, this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   void ChandrasekharOnsetSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<EnergyFile<LinearThermalConvectionTraits::CodType> > pCodEnergy(new EnergyFile<LinearThermalConvectionTraits::CodType>(this->mCodC, "cod", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<LinearThermalConvectionTraits::VelType> > pVelEnergy(new EnergyFile<LinearThermalConvectionTraits::VelType>(this->mVelV, "vel", this->mSimControl.tsParams(),1));

      // Add codensity energy to ASCII output
      this->mIOSys.addASCIIWriter(pCodEnergy);
      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<SpectrumFile<LinearThermalConvectionTraits::CodType> > pCodSpectrum(new SpectrumFile<LinearThermalConvectionTraits::CodType>(this->mCodC, "cod"));
      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<LinearThermalConvectionTraits::VelType> > pVelSpectrum(new SpectrumFile<LinearThermalConvectionTraits::VelType>(this->mVelV, "vel"));

      // Add codensity energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pCodSpectrum);
      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}
