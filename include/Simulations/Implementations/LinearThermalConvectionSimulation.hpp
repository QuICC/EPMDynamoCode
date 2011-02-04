/** \file LinearThermalConvectionSimulation.hpp
 *  \brief Implementation of a linear thermal convection simulation
 */

#ifndef LINEARTHERMALCONVECTIONSIMULATION_HPP
#define LINEARTHERMALCONVECTIONSIMULATION_HPP

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
#include "Simulations/Traits/LinearThermalConvectionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

#include "Equations/Transport/TransportLinear.hpp"
#include "Equations/NavierStokes/NavierStokesLinThermal.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a linear thermal convection simulation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class LinearThermalConvectionSimulation: public SimulationBase<TSimType>
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~LinearThermalConvectionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         LinearThermalConvectionSimulation();

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
          * @brief Codensity scalar
          */
         typename LinearThermalConvectionTraits<TSimType>::CodType   mCodC;

         /**
          * @brief Velocity field
          */
         typename LinearThermalConvectionTraits<TSimType>::VelType   mVelV;

         /**
          * @brief Transport equation
          */
         TransportLinear<TSimType, LinearThermalConvectionTraits>    mTransport;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesLinThermal<TSimType, LinearThermalConvectionTraits>    mNavierStokes;
   };

   template <typename TSimType> LinearThermalConvectionSimulation<TSimType>::LinearThermalConvectionSimulation()
      : mCodC(this->mpTrunc, this->mTransform), mVelV(this->mpTrunc, this->mTransform), mTransport(mCodC, mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mNavierStokes(mVelV, mCodC, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::initEquations()
   {
      // Set boundary condition to the transport equation
      SmartBC  pZeroBC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
      this->mTransport.addBC(pZeroBC);

      // Set boundary condition to the Navier-Stokes equation
      SmartBC  pSFBC(new StressFreeTorBC<TSimType>(this->mTransform.radBasis()));
      SmartBC  pDDBC(new DDRadialBC<TSimType>(this->mTransform.radBasis()));
      this->mNavierStokes.addTorBC(pSFBC);
      // Order of Poloidal BCs is important
      this->mNavierStokes.addPolBC(pZeroBC);
      this->mNavierStokes.addPolBC(pDDBC);

      // Initialise the transport equation
      this->mTransport.init();

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();

      TimestepConfig::overruleMaxTimestep(1e-5);
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::configureTransforms()
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

      // Configure the SSH manipulator
      this->configureSHManipulator();

      //
      // Configure transform nesting
      //

      // Configure the transforms' nesting setup
      this->configureTransformNesting();
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::updateEquationsRTP(const int step)
   {
      // Update RTP values of the transport equation
      this->mTransport.updateRTP(step);

      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::updateEquationsRHS()
   {
      // Update RHS of the transport equation
      this->mTransport.updateRHS();

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateRTPCFLTimestep(this->mVelV.oc().rtp());
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::transformEquationsRHS(const int step)
   {
      // Update RHS of the transport equation
      this->mTransport.transformRHS(step);

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::addExternalInfluence()
   {
      // Add Thermal source term
      this->mTransport.addSourceTerm();
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::timestepEquations()
   {
      // Timestep the transport equation
      this->mTransport.timestep();

      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<TSimType, LinearThermalConvectionTraits> > pInState(new StateFileReader<TSimType, LinearThermalConvectionTraits>(this->mCodC, this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create a source file reader for the codensity source
      EPMSHARED_PTR<SourceFileReader<TSimType, LinearThermalConvectionTraits> > pSource(new SourceFileReader<TSimType, LinearThermalConvectionTraits>(this->mCodC));

      // Read in source state
      this->mIOSys.useSource(pSource);

      // Rescale codensity source to real value from parameters
      this->mCodC.rOcSrc().rescale(this->mEqParams.codSourceScale());
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<TSimType, LinearThermalConvectionTraits> >  pOutState(new StateFileWriter<TSimType, LinearThermalConvectionTraits>(this->mCodC, this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   template <typename TSimType> void LinearThermalConvectionSimulation<TSimType>::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<EnergyFile<TSimType, typename LinearThermalConvectionTraits<TSimType>::CodType> > pCodEnergy(new EnergyFile<TSimType, typename LinearThermalConvectionTraits<TSimType>::CodType>(this->mCodC, "cod", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<TSimType, typename LinearThermalConvectionTraits<TSimType>::VelType> > pVelEnergy(new EnergyFile<TSimType, typename LinearThermalConvectionTraits<TSimType>::VelType>(this->mVelV, "vel", this->mSimControl.tsParams()));

      // Add codensity energy to ASCII output
      this->mIOSys.addASCIIWriter(pCodEnergy);
      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<SpectrumFile<TSimType, typename LinearThermalConvectionTraits<TSimType>::CodType> > pCodSpectrum(new SpectrumFile<TSimType, typename LinearThermalConvectionTraits<TSimType>::CodType>(this->mCodC, "cod"));
      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<TSimType, typename LinearThermalConvectionTraits<TSimType>::VelType> > pVelSpectrum(new SpectrumFile<TSimType, typename LinearThermalConvectionTraits<TSimType>::VelType>(this->mVelV, "vel"));

      // Add codensity energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pCodSpectrum);
      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}

#endif // LINEARTHERMALCONVECTIONSIMULATION_HPP
