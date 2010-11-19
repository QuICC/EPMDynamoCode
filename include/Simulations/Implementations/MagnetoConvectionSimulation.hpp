/** \file MagnetoConvectionSimulation.hpp
 *  \brief Implementation of the magneto-convection simulation
 */

#ifndef MAGNETOCONVECTIONSIMULATION_HPP
#define MAGNETOCONVECTIONSIMULATION_HPP

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
#include "Simulations/Traits/MagnetoConvectionTraits.hpp"

#include "IO/HDF5/Imposed/ImposedFieldReader.hpp"
#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"

#include "Equations/Induction/InductionMHD.hpp"
#include "Equations/NavierStokes/NavierStokesMagnetic.hpp"
#include "Equations/Transport/TransportMHD.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/InsulatingBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the magneto-convection simulation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class MagnetoConvectionSimulation: public SimulationBase<TSimType>
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~MagnetoConvectionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         MagnetoConvectionSimulation();

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
         typename MagnetoConvectionTraits<TSimType>::CodType   mCodC;

         /**
          * @brief Magnetic field
          */
         typename MagnetoConvectionTraits<TSimType>::MagType   mMagB;

         /**
          * @brief Velocity field
          */
         typename MagnetoConvectionTraits<TSimType>::VelType   mVelV;

         /**
          * @brief Induction equation
          */
         InductionMHD<TSimType, MagnetoConvectionTraits>    mInduction;

         /**
          * @brief Transport equation
          */
         TransportMHD<TSimType, MagnetoConvectionTraits>    mTransport;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesMagnetic<TSimType, MagnetoConvectionTraits>    mNavierStokes;
   };

   template <typename TSimType> MagnetoConvectionSimulation<TSimType>::MagnetoConvectionSimulation()
      : mCFLFactor(1.0), mCodC(this->mpTrunc, this->mTransform), mMagB(this->mpTrunc, this->mTransform), mVelV(this->mpTrunc, this->mTransform), mInduction(mMagB, mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mTransport(mCodC, mVelV, this->mTransform, this->mSimControl.tsParams(), this->mEqParams), mNavierStokes(mVelV, mMagB, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
      // Set the CFL factor to L*(L+1)
      int l = this->mpTrunc->sim()->hoz()->nL();
      this->mCFLFactor = static_cast<EPMFloat>(l*(l+1));
   }

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::initEquations()
   {
      // Set boundary condition to the transport equation
      SmartBC  pZeroBCC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
      //this->mTransport.addBC(pDiffHeatBC);
      this->mTransport.addBC(pZeroBCC);

      // Set boundary condition to the Induction equation
      SmartBC  pZeroBC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
      SmartBC  pInsulatingBC(new InsulatingBC<TSimType>(this->mTransform.radBasis()));
      mInduction.addTorBC(pZeroBC);
      mInduction.addPolBC(pInsulatingBC);

      // Set boundary condition to the Navier-Stokes equation
      //SmartBC  pSFBC(new StressFreeTorBC<TSimType>(this->mTransform.radBasis()));
      SmartBC  pNSBC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
      //SmartBC  pDDBC(new DDRadialBC<TSimType>(this->mTransform.radBasis()));
      SmartBC  pDBC(new DRadialBC<TSimType>(this->mTransform.radBasis()));
      this->mNavierStokes.addTorBC(pNSBC);
      // Order of Poloidal BCs is important
      this->mNavierStokes.addPolBC(pZeroBC);
      this->mNavierStokes.addPolBC(pDBC);

      // Initialise the induction equation
      this->mInduction.init();

      // Initialise the transport equation
      this->mTransport.init();

      // Initialise the Navier-Stokes equation
      this->mNavierStokes.init();
   }

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::updateEquationsRTP(const int step)
   {
      // Update RTP values of the induction equation
      this->mInduction.updateRTP(step);

      // Update RTP values of the transport equation
      this->mTransport.updateRTP(step);

      // Update RTP values of the Navier-Stokes equation
      this->mNavierStokes.updateRTP(step);
   }

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::updateEquationsRHS()
   {
      // Update RHS of the induction equation
      this->mInduction.updateRHS();

      // Update RHS of the transport equation
      this->mTransport.updateRHS();

      // Update RHS of the Navier-Stokes equation
      this->mNavierStokes.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateCFLTimestep(mMagB.oc().rtp(), mVelV.oc().rtp(), mCFLFactor);
   }

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::transformEquationsRHS(const int step)
   {
      // Transform RHS of the induction equation
      mInduction.transformRHS(step);

      // Transform RHS of the transport equation
      this->mTransport.transformRHS(step);

      // Transform RHS of the Navier-Stokes equation
      mNavierStokes.transformRHS(step);
   }

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::addExternalInfluence()
   {
      // Add Thermal source term
      this->mTransport.addSourceTerm();
   }

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::timestepEquations()
   {
      // Timestep the induction equation
      mInduction.timestep();

      // Timestep the transport equation
      this->mTransport.timestep();

      // Timestep the Navier-Stokes equation
      mNavierStokes.timestep();
   }

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::configureTransforms()
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

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::initFields()
   {
      // Read in initial state
      EPMSHARED_PTR<StateFileReader<TSimType, MagnetoConvectionTraits> > pInState(new StateFileReader<TSimType, MagnetoConvectionTraits>(this->mCodC, this->mMagB, this->mVelV,  "_initial"));

      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create a source file reader for the codensity source
      EPMSHARED_PTR<SourceFileReader<TSimType, MagnetoConvectionTraits> > pSource(new SourceFileReader<TSimType, MagnetoConvectionTraits>(this->mCodC));

      // Read in source state
      this->mIOSys.useSource(pSource);

      // Rescale codensity source
      this->mCodC.rOcSrc().rescale(this->mEqParams.codSourceScale());

      // Read in imposed magnetic field
      EPMSHARED_PTR<ImposedFieldReader<TSimType, MagnetoConvectionTraits> > pIField(new ImposedFieldReader<TSimType, MagnetoConvectionTraits>(this->mMagB));

      this->mIOSys.useImposedField(pIField);

      // Scale imposed field to correct values
      this->mMagB.rOc().rImposed().rescale(this->mEqParams.imposedMagScale());

      // Set the energy scale for the magnetic field
      this->mMagB.rOc().rPerturbation().setEnergyScale(this->mEqParams.meFactor());

      // Set the energy scale for the velocity field
      this->mVelV.rOc().rPerturbation().setEnergyScale(this->mEqParams.keFactor());
   }

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<TSimType, MagnetoConvectionTraits> >  pOutState(new StateFileWriter<TSimType, MagnetoConvectionTraits>(this->mCodC, this->mMagB, this->mVelV, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   template <typename TSimType> void MagnetoConvectionSimulation<TSimType>::addASCIIOutput()
   {
      // Create a energy ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<EnergyFile<TSimType, typename MagnetoConvectionTraits<TSimType>::CodType> > pCodEnergy(new EnergyFile<TSimType, typename MagnetoConvectionTraits<TSimType>::CodType>(this->mCodC, "cod", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<EnergyFile<TSimType, typename MagnetoConvectionTraits<TSimType>::MagType> > pMagEnergy(new EnergyFile<TSimType, typename MagnetoConvectionTraits<TSimType>::MagType>(this->mMagB, "mag", this->mSimControl.tsParams()));
      // Create a energy ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<EnergyFile<TSimType, typename MagnetoConvectionTraits<TSimType>::VelType> > pVelEnergy(new EnergyFile<TSimType, typename MagnetoConvectionTraits<TSimType>::VelType>(this->mVelV, "vel", this->mSimControl.tsParams()));

      // Add codensity energy to ASCII output
      this->mIOSys.addASCIIWriter(pCodEnergy);
      // Add magnetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pMagEnergy);
      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pVelEnergy);

      // Create a energy spectrum ASCII diagnostic file for the codensity scalar
      EPMSHARED_PTR<SpectrumFile<TSimType, typename MagnetoConvectionTraits<TSimType>::CodType> > pCodSpectrum(new SpectrumFile<TSimType, typename MagnetoConvectionTraits<TSimType>::CodType>(this->mCodC, "cod"));
      // Create a energy spectrum ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<SpectrumFile<TSimType, typename MagnetoConvectionTraits<TSimType>::MagType> > pMagSpectrum(new SpectrumFile<TSimType, typename MagnetoConvectionTraits<TSimType>::MagType>(this->mMagB, "mag"));
      // Create a energy spectrum ASCII diagnostic file for the velocity field
      EPMSHARED_PTR<SpectrumFile<TSimType, typename MagnetoConvectionTraits<TSimType>::VelType> > pVelSpectrum(new SpectrumFile<TSimType, typename MagnetoConvectionTraits<TSimType>::VelType>(this->mVelV, "vel"));

      // Add codensity energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pCodSpectrum);
      // Add magnetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pMagSpectrum);
      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pVelSpectrum);
   }

}

#endif // MAGNETOCONVECTIONSIMULATION_HPP
