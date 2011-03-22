/** \file VelocityRotDiffusionSimulation.hpp
 *  \brief Implementation of a velocity diffusion simulation
 */

#ifndef VELOCITYROTDIFFUSIONSIMULATION_HPP
#define VELOCITYROTDIFFUSIONSIMULATION_HPP

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
#include "Simulations/Traits/VelocityRotDiffusionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

#include "Equations/NavierStokes/NavierStokesRotDiffusion.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a velocity diffusion simulation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class VelocityRotDiffusionSimulation: public SimulationBase<TSimType>
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~VelocityRotDiffusionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         VelocityRotDiffusionSimulation();

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
         typename VelocityRotDiffusionTraits<TSimType>::VelType   mVelV;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesRotDiffusion<TSimType, VelocityRotDiffusionTraits>    mNavierStokes;
   };

   template <typename TSimType> VelocityRotDiffusionSimulation<TSimType>::VelocityRotDiffusionSimulation()
      : mVelV(this->mpTrunc, this->mTransform), mNavierStokes(mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::initEquations()
   {
      // Create Zero boundary condition pointer
      SmartBC  pZeroBC(new ZeroBC<TSimType>(this->mTransform.radBasis()));

      // Set boundary condition to the Navier-Stokes equation
      if(this->mIOSys.cfg()->aBC()(1) == 1)
      {
         SmartBC  pSFBC(new StressFreeTorBC<TSimType>(this->mTransform.radBasis()));
         SmartBC  pDDBC(new DDRadialBC<TSimType>(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pSFBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDDBC);
      } else
      {
         SmartBC  pNSBC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
         SmartBC  pDBC(new DRadialBC<TSimType>(this->mTransform.radBasis()));

         // Toroidal velocity BC
         this->mNavierStokes.addTorBC(pNSBC);

         // Order of Poloidal BCs is important
         this->mNavierStokes.addPolBC(pZeroBC);
         this->mNavierStokes.addPolBC(pDBC);
      }

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::configureTransforms()
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

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::updateEquationsRTP(const int step)
   {
      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::updateEquationsRHS()
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mVelV.oc().rtp());
   }

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::transformEquationsRHS(const int step)
   {
      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::addExternalInfluence()
   {
   }

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::timestepEquations()
   {
      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<TSimType, VelocityRotDiffusionTraits> > pInState(new StateFileReader<TSimType, VelocityRotDiffusionTraits>(this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());
   }

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<TSimType, VelocityRotDiffusionTraits> >  pOutState(new StateFileWriter<TSimType, VelocityRotDiffusionTraits>(this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   template <typename TSimType> void VelocityRotDiffusionSimulation<TSimType>::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<TSimType, typename VelocityRotDiffusionTraits<TSimType>::VelType> > pVelEnergy(new EnergyFile<TSimType, typename VelocityRotDiffusionTraits<TSimType>::VelType>(mVelV, "vel", this->mSimControl.tsParams()));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<TSimType, typename VelocityRotDiffusionTraits<TSimType>::VelType> > pVelSpectrum(new SpectrumFile<TSimType, typename VelocityRotDiffusionTraits<TSimType>::VelType>(this->mVelV, "vel"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}

#endif // VELOCITYROTDIFFUSIONSIMULATION_HPP
