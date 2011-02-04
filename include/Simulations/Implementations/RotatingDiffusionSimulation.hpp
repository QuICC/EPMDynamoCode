/** \file RotatingDiffusionSimulation.hpp
 *  \brief Implementation of a rotating velocity diffusion simulation
 */

#ifndef ROTATINGDIFFUSIONSIMULATION_HPP
#define ROTATINGDIFFUSIONSIMULATION_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Simulations/SimulationBase.hpp"
#include "Simulations/Traits/RotatingDiffusionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

#include "Equations/NavierStokes/NavierStokesDiffusion.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a rotating velocity diffusion simulation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class RotatingDiffusionSimulation: public SimulationBase<TSimType>
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~RotatingDiffusionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         RotatingDiffusionSimulation();

         /**
          * @brief Initialise the fields
          */
         void initFields();

         /**
          * @brief Add HDF5 output files
          */
         void addHDF5Output();

         /**
          * @brief Add ASCII output files
          */
         void addASCIIOutput();

         /**
          * @brief Initialise the equations
          */
         void initEquations();

         /**
          * @brief Update equations RTP values
          *
          * \param step Current step in a multistep transform
          */
         void updateEquationsRTP(const int step);

         /**
          * @brief Update equations RHS
          */
         void updateEquationsRHS();

         /**
          * @brief Transform equations RHS
          *
          * \param step Current step in a multistep transform
          */
         void transformEquationsRHS(const int step);

         /**
          * @brief Include external influence
          */
         void addExternalInfluence();

         /**
          * @brief Timestep the equations
          */
         void timestepEquations();

         /**
          * @brief Configure the transforms
          */
         void configureTransforms();

      private:
         /**
          * @brief Velocity field
          */
         typename RotatingDiffusionTraits<TSimType>::VelType   mVelV;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesDiffusion<TSimType, RotatingDiffusionTraits>    mNavierStokes;
   };

   template <typename TSimType> RotatingDiffusionSimulation<TSimType>::RotatingDiffusionSimulation()
      : mVelV(this->mpTrunc, this->mTransform), mNavierStokes(mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::initEquations()
   {
      // Create Zero boundary condition pointer
      SmartBC  pZeroBC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
      // Create toroidal stress free boundary condition pointer
      SmartBC  pSFBC(new StressFreeTorBC<TSimType>(this->mTransform.radBasis()));
      // Create poloidal stress free boundary condition pointer
      SmartBC  pDDBC(new DDRadialBC<TSimType>(this->mTransform.radBasis()));

      // Set boundary condition for the toroidal component
      this->mNavierStokes.addTorBC(pSFBC);

      // Set boundary conditions for the Poloidal component
      // Order of BCs is IMPORTANT (first one is applied to influence matrix steps)
      this->mNavierStokes.addPolBC(pZeroBC);
      this->mNavierStokes.addPolBC(pDDBC);

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::configureTransforms()
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

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::updateEquationsRTP(const int step)
   {
      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::updateEquationsRHS()
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mVelV.oc().rtp());
   }

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::transformEquationsRHS(const int step)
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::addExternalInfluence()
   {
   }

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::timestepEquations()
   {
      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<TSimType, RotatingDiffusionTraits> > pInState(new StateFileReader<TSimType, RotatingDiffusionTraits>(this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());
   }

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<TSimType, RotatingDiffusionTraits> >  pOutState(new StateFileWriter<TSimType, RotatingDiffusionTraits>(this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   template <typename TSimType> void RotatingDiffusionSimulation<TSimType>::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<TSimType, typename RotatingDiffusionTraits<TSimType>::VelType> > pVelEnergy(new EnergyFile<TSimType, typename RotatingDiffusionTraits<TSimType>::VelType>(mVelV, "vel", this->mSimControl.tsParams()));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<TSimType, typename RotatingDiffusionTraits<TSimType>::VelType> > pVelSpectrum(new SpectrumFile<TSimType, typename RotatingDiffusionTraits<TSimType>::VelType>(this->mVelV, "vel"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}

#endif // ROTATINGDIFFUSIONSIMULATION_HPP
