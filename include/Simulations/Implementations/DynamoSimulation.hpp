/** \file DynamoSimulation.hpp
 *  \brief Implementation of a dynamo simulation
 */

#ifndef DYNAMOSIMULATION_HPP
#define DYNAMOSIMULATION_HPP

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
#include "Simulations/Traits/DynamoTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"

#include "Equations/Induction/InductionMHD.hpp"
#include "Equations/Transport/TransportMHD.hpp"
#include "Equations/NavierStokes/NavierStokesMHD.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/InsulatingBC.hpp"
#include "BoundaryConditions/L0Harmonic/DifferentialHeatingBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a dynamo simulation as a simulation policy class
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class DynamoSimulation: public SimulationBase<TSimType>
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~DynamoSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         DynamoSimulation();

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
          * @brief Angular distance factor for CFL condition
          */
         EPMFloat mCFLFactor;

         /**
          * @brief Codensity scalar
          */
         typename DynamoTraits<TSimType>::CodType   mCodC;

         /**
          * @brief Magnetic field
          */
         typename DynamoTraits<TSimType>::MagType   mMagB;

         /**
          * @brief Velocity field
          */
         typename DynamoTraits<TSimType>::VelType   mVelV;

         /**
          * @brief Induction equation
          */
         InductionMHD<TSimType, DynamoTraits>    mInduction;

         /**
          * @brief Transport equation
          */
         TransportMHD<TSimType, DynamoTraits>    mTransport;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesMHD<TSimType, DynamoTraits>    mNavierStokes;
   };

   template <typename TSimType> DynamoSimulation<TSimType>::DynamoSimulation()
      : mCFLFactor(1.0), mCodC(this->mpTrunc, this->mTransform), mMagB(this->mpTrunc, this->mTransform), mVelV(this->mpTrunc, this->mTransform), mInduction(mMagB, mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mTransport(mCodC, mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mNavierStokes(mVelV, mMagB, mCodC, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
      // Set the CFL factor to L*(L+1)
      int l = this->mpTrunc->sim()->hoz()->nL();
      this->mCFLFactor = static_cast<EPMFloat>(l*(l+1));
   }

   template <typename TSimType> void DynamoSimulation<TSimType>::initEquations()
   {
      // Set boundary condition to the transport equation
      SmartBC  pZeroBCC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
      //this->mTransport.addBC(pDiffHeatBC);
      this->mTransport.addBC(pZeroBCC);

      // Set boundary condition to the Induction equation
      SmartBC  pZeroBC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
      SmartBC  pInsulatingBC(new InsulatingBC<TSimType>(this->mTransform.radBasis()));
      this->mInduction.addTorBC(pZeroBC);
      this->mInduction.addPolBC(pInsulatingBC);

      // Set boundary condition to the Navier-Stokes equation
      SmartBC  pSFBC(new StressFreeTorBC<TSimType>(this->mTransform.radBasis()));
      SmartBC  pDDBC(new DDRadialBC<TSimType>(this->mTransform.radBasis()));
      this->mNavierStokes.addTorBC(pSFBC);
      // Order of Poloidal BCs is important
      this->mNavierStokes.addPolBC(pZeroBC);
      this->mNavierStokes.addPolBC(pDDBC);

      // Initialise the induction equation
      this->mInduction.init();

      // Initialise the transport equation
      this->mTransport.init();

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   template <typename TSimType> void DynamoSimulation<TSimType>::updateEquationsRTP(const int step)
   {
      // Update RTP values of the induction equation
      this->mInduction.updateRTP(step);

      // Update RTP values of the transport equation
      this->mTransport.updateRTP(step);

      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   template <typename TSimType> void DynamoSimulation<TSimType>::updateEquationsRHS()
   {
      // Update RHS of the induction equation
      this->mInduction.updateRHS();

      // Update RHS of the transport equation
      this->mTransport.updateRHS();

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateCFLTimestep(this->mMagB.oc().rtp(), this->mVelV.oc().rtp(), this->mCFLFactor);
   }

   template <typename TSimType> void DynamoSimulation<TSimType>::transformEquationsRHS(const int step)
   {
      // Transform RHS of the induction equation
      this->mInduction.transformRHS(step);

      // Transform RHS of the transport equation
      this->mTransport.transformRHS(step);

      // Transform RHS of the Navier-Stokes equation
      this->mNavierStokes.transformRHS(step);
   }

   template <typename TSimType> void DynamoSimulation<TSimType>::addExternalInfluence()
   {
      // Add Thermal source term
      this->mTransport.addSourceTerm();
   }

   template <typename TSimType> void DynamoSimulation<TSimType>::timestepEquations()
   {
      // Timestep the induction equation
      this->mInduction.timestep();

      // Timestep the transport equation
      this->mTransport.timestep();

      // Timestep the Navier-Stokes equation
      this->mNavierStokes.timestep();
   }

   template <typename TSimType> void DynamoSimulation<TSimType>::configureTransforms()
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

   template <typename TSimType> void DynamoSimulation<TSimType>::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<TSimType, DynamoTraits> > pInState(new StateFileReader<TSimType, DynamoTraits>(this->mCodC, this->mMagB, this->mVelV,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create a source file reader for the codensity source
      EPMSHARED_PTR<SourceFileReader<TSimType, DynamoTraits> > pSource(new SourceFileReader<TSimType, DynamoTraits>(this->mCodC));

      // Read in source state
      this->mIOSys.useSource(pSource);
   }

   template <typename TSimType> void DynamoSimulation<TSimType>::addHDF5Output()
   {
      // Create a state file writer for regular state output
      EPMSHARED_PTR<StateFileWriter<TSimType, DynamoTraits> >  pOutState(new StateFileWriter<TSimType, DynamoTraits>(this->mCodC, this->mMagB, this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      // Add state file writer to hdf5 outputs
      this->mIOSys.addHDF5Writer(pOutState);
   }

   template <typename TSimType> void DynamoSimulation<TSimType>::addASCIIOutput()
   {
      // Create a energy ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<EnergyFile<TSimType, typename DynamoTraits<TSimType>::CodType> > pCodEnergy(new EnergyFile<TSimType, typename DynamoTraits<TSimType>::CodType>(this->mCodC, "cod", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<EnergyFile<TSimType, typename DynamoTraits<TSimType>::MagType> > pMagEnergy(new EnergyFile<TSimType, typename DynamoTraits<TSimType>::MagType>(this->mMagB, "mag", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<TSimType, typename DynamoTraits<TSimType>::VelType> > pVelEnergy(new EnergyFile<TSimType, typename DynamoTraits<TSimType>::VelType>(this->mVelV, "vel", this->mSimControl.tsParams()));

      // Add codensity energy to ASCII output
      this->mIOSys.addASCIIWriter(pCodEnergy);
      // Add magnetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pMagEnergy);
      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<SpectrumFile<TSimType, typename DynamoTraits<TSimType>::CodType> > pCodSpectrum(new SpectrumFile<TSimType, typename DynamoTraits<TSimType>::CodType>(this->mCodC, "cod"));
      // Create a energy spectrum ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<SpectrumFile<TSimType, typename DynamoTraits<TSimType>::MagType> > pMagSpectrum(new SpectrumFile<TSimType, typename DynamoTraits<TSimType>::MagType>(this->mMagB, "mag"));
      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<TSimType, typename DynamoTraits<TSimType>::VelType> > pVelSpectrum(new SpectrumFile<TSimType, typename DynamoTraits<TSimType>::VelType>(this->mVelV, "vel"));

      // Add codensity energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pCodSpectrum);
      // Add magnetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pMagSpectrum);
      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}

#endif // DYNAMOSIMULATION_HPP
