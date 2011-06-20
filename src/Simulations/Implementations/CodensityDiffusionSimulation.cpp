/** \file CodensityDiffusionSimulation.cpp
 *  \brief Implementation of a codensity diffusion simulation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Simulations/Implementations/CodensityDiffusionSimulation.hpp"

// Project includes
//
#include "BoundaryConditions/L0Harmonic/HeatFluxBC.hpp"

namespace EPMDynamo {

   CodensityDiffusionSimulation::CodensityDiffusionSimulation()
      : mCodC(this->mpTrunc, this->mTransform), mTransport(this->mCodC, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void CodensityDiffusionSimulation::initEquations()
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

      // Initialise the transport equation
      this->mTransport.init();
   }

   void CodensityDiffusionSimulation::configureTransforms()
   {
      //
      // Setup the SSH transform data manipulator
      //

      // Register transform data packs of transport equation for SSH
      this->registerSSHPacks(this->mTransport.nFSSHPacks(), this->mTransport.nSSHBPacks());

      // Configure the SSH manipulator
      this->configureSSHManipulator();

      //
      // Setup the SH transform data manipulator
      //

      // Register transform data packs of transport equation for SH
      this->registerSHPacks(this->mTransport.nFSHPacks(), this->mTransport.nSHBPacks());

      // Configure the SSH manipulator
      this->configureSHManipulator();

      //
      // Configure transform nesting
      //

      // Configure the transforms' nesting setup
      this->configureTransformNesting();
   }

   void CodensityDiffusionSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the transport equation
      this->mTransport.updateRTP(step);
   }

   void CodensityDiffusionSimulation::updateEquationsRHS()
   {
      // Update RHS of the transport equation
      this->mTransport.updateRHS();
   }

   void CodensityDiffusionSimulation::transformEquationsRHS(const int step)
   {
      // Update RHS of the transport equation
      this->mTransport.transformRHS(step);
   }

   void CodensityDiffusionSimulation::addExternalInfluence()
   {
   }

   void CodensityDiffusionSimulation::timestepEquations()
   {
      // Timestep the transport equation
      this->mTransport.timestep();
   }

   void CodensityDiffusionSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<CodensityDiffusionTraits> > pInState(new StateFileReader<CodensityDiffusionTraits>(this->mCodC,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create a source file reader for the codensity source
      EPMSHARED_PTR<SourceFileReader<CodensityDiffusionTraits> > pSource(new SourceFileReader<CodensityDiffusionTraits>(this->mCodC));

      // Read in source state
      this->mIOSys.useSource(pSource);

      // Rescale the source term to correct parameter related value
      this->mCodC.rOcSrc().rescale(this->mEqParams.codSourceScale());
   }

   void CodensityDiffusionSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<CodensityDiffusionTraits> >  pOutState(new StateFileWriter<CodensityDiffusionTraits>(this->mCodC, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   void CodensityDiffusionSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the codensity field
      EPMSHARED_PTR<EnergyFile<CodensityDiffusionTraits::CodType> > pCodEnergy(new EnergyFile<CodensityDiffusionTraits::CodType>(mCodC, "cod", this->mSimControl.tsParams()));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pCodEnergy);

      // Create a energy spectrum ASCII diagnostic file for the codensity field
      EPMSHARED_PTR<SpectrumFile<CodensityDiffusionTraits::CodType> > pCodSpectrum(new SpectrumFile<CodensityDiffusionTraits::CodType>(this->mCodC, "cod"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pCodSpectrum);
   }

}
