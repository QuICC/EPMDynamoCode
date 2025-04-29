/** \file KinematicInductionSimulation.cpp
 *  \brief Implementation of a kinematic induction simulation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Simulations/Implementations/KinematicInductionSimulation.hpp"

// Project includes
//

namespace EPMDynamo {

   KinematicInductionSimulation::KinematicInductionSimulation()
      : mMagB(this->mpTrunc, this->mTransform), mVelV(this->mpTrunc, this->mTransform), mInduction(this->mMagB, this->mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   void KinematicInductionSimulation::initEquations()
   {
      // Set boundary condition to the induction equation
      if(this->mIOSys.cfg()->aBC()(2) == 0)
      {
         // Set insulator boundary condition
         SmartBC  pZeroBC(new ZeroBC(this->mTransform.radBasis()));
         SmartBC  pInsulatingBC(new InsulatingBC(this->mTransform.radBasis()));

         this->mInduction.addTorBC(pZeroBC);
         this->mInduction.addPolBC(pInsulatingBC);
      } else if(this->mIOSys.cfg()->aBC()(2) == 1)
      {
         // Set conductor boundary condition
         SmartBC  pConductorTorBC(new ConductorTorBC(this->mTransform.radBasis()));
         SmartBC  pConductorPolBC(new ConductorPolBC(this->mTransform.radBasis()));

         this->mInduction.addTorBC(pConductorTorBC);
         this->mInduction.addPolBC(pConductorPolBC);
      } else
      {
         throw EPMException("KinematicInductionSimulation::initEquations", "Did not know what to do with Magnetic BC");
      }

      // Initialise the Induction equation
      this->mInduction.init();
   }

   void KinematicInductionSimulation::configureTransforms()
   {
      //
      // Setup the SSH transform data manipulator
      //

      // Register transform data packs of induction equation for SSH
      this->registerSSHPacks(this->mInduction.nFSSHPacks(), this->mInduction.nSSHBPacks());

      // Configure the SSH manipulator
      this->configureSSHManipulator();

      //
      // Setup the SH transform data manipulator
      //

      // Register transform data packs of induction equation for SH
      this->registerSHPacks(this->mInduction.nFSHPacks(), this->mInduction.nSHBPacks());

      // Configure the SSH manipulator
      this->configureSHManipulator();

      //
      // Configure transform nesting
      //

      // Configure the transforms' nesting setup
      this->configureTransformNesting();
   }

   void KinematicInductionSimulation::updateEquationsRTP(const int step)
   {
      // Update RTP values of the Induction equation
      this->mInduction.updateRTP(step);
   }

   void KinematicInductionSimulation::updateEquationsRHS()
   {
      // Update RHS of the induction equation
      this->mInduction.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mMagB.oc().rtp(), this->mVelV.oc().rtp());
   }

   void KinematicInductionSimulation::transformEquationsRHS(const int step)
   {
      // Update RHS of the induction equation
      this->mInduction.transformRHS(step);
   }

   void KinematicInductionSimulation::addExternalInfluence()
   {
   }

   void KinematicInductionSimulation::timestepEquations()
   {
      // Timestep the induction equation
      this->mInduction.timestep();
   }

   void KinematicInductionSimulation::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<KinematicInductionTraits> > pInState(new StateFileReader<KinematicInductionTraits>(this->mMagB, this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Set the energy scale for the magnetic field
      this->mMagB.rOc().rPerturbation().setEnergyScale(this->mEqParams.meFactor());

      // Set the energy scale for the magnetic field
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());

      // Compute physical values for velocity field
      for(int i=0; i < this->mTransformSteps; ++i)
      {
         this->combineRTPTransforms(i);

         this->mVelV.rOc().transform(i);
      }
   }

   void KinematicInductionSimulation::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<KinematicInductionTraits> >  pOutState(new StateFileWriter<KinematicInductionTraits>(this->mMagB, this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   void KinematicInductionSimulation::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<EnergyFile<KinematicInductionTraits::MagType> > pMagEnergy(new EnergyFile<KinematicInductionTraits::MagType>(this->mMagB, "mag", this->mSimControl.tsParams(),1));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pMagEnergy);

      // Create a energy spectrum ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<SpectrumFile<KinematicInductionTraits::MagType> > pMagSpectrum(new SpectrumFile<KinematicInductionTraits::MagType>(this->mMagB, "mag"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pMagSpectrum);
   }

}
