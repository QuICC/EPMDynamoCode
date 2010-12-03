/** \file MagneticDiffusionSimulation.hpp
 *  \brief Implementation of a magnetic diffusion simulation
 */

#ifndef MAGNETICDIFFUSIONSIMULATION_HPP
#define MAGNETICDIFFUSIONSIMULATION_HPP

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
#include "Simulations/Traits/MagneticDiffusionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

#include "Equations/Induction/InductionDiffusion.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/InsulatingBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a magnetic diffusion simulation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class MagneticDiffusionSimulation: public SimulationBase<TSimType>
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~MagneticDiffusionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         MagneticDiffusionSimulation();

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
          * @brief Magnetic field
          */
         typename MagneticDiffusionTraits<TSimType>::MagType   mMagB;

         /**
          * @brief Induction equation
          */
         InductionDiffusion<TSimType, MagneticDiffusionTraits>    mInduction;
   };

   template <typename TSimType> MagneticDiffusionSimulation<TSimType>::MagneticDiffusionSimulation()
      : mCFLFactor(1.0), mMagB(this->mpTrunc, this->mTransform), mInduction(mMagB, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
      // Set the CFL factor to L*(L+1)
      int l = this->mpTrunc->sim()->hoz()->nL();
      this->mCFLFactor = static_cast<EPMFloat>(l*(l+1));
   }

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::initEquations()
   {
      // Set boundary condition to the Induction equation
      SmartBC  pZeroBC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
      SmartBC  pInsulatingBC(new InsulatingBC<TSimType>(this->mTransform.radBasis()));
      this->mInduction.addTorBC(pZeroBC);
      this->mInduction.addPolBC(pInsulatingBC);

      // Initialise the Induction equation
      this->mInduction.init();
   }

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::configureTransforms()
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

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::updateEquationsRTP(const int step)
   {
      // Update RTP values of the Induction equation
      this->mInduction.updateRTP(step);
   }

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::updateEquationsRHS()
   {
      // Update RHS of the induction equation
      this->mInduction.updateRHS();

      // Update the CFL timestep condition
      this->mSimControl.tsControl().updateCFLTimestep(this->mMagB.oc().rtp(), this->mCFLFactor);
   }

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::transformEquationsRHS(const int step)
   {
      // Update RHS of the induction equation
      this->mInduction.transformRHS(step);
   }

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::addExternalInfluence()
   {
   }

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::timestepEquations()
   {
      // Timestep the induction equation
      this->mInduction.timestep();
   }

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<TSimType, MagneticDiffusionTraits> > pInState(new StateFileReader<TSimType, MagneticDiffusionTraits>(this->mMagB,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());
   }

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<TSimType, MagneticDiffusionTraits> >  pOutState(new StateFileWriter<TSimType, MagneticDiffusionTraits>(this->mMagB, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   template <typename TSimType> void MagneticDiffusionSimulation<TSimType>::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<EnergyFile<TSimType, typename MagneticDiffusionTraits<TSimType>::MagType> > pMagEnergy(new EnergyFile<TSimType, typename MagneticDiffusionTraits<TSimType>::MagType>(mMagB, "mag", this->mSimControl.tsParams()));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pMagEnergy);

      // Create a energy spectrum ASCII diagnostic file for the magnetic field
      EPMSHARED_PTR<SpectrumFile<TSimType, typename MagneticDiffusionTraits<TSimType>::MagType> > pMagSpectrum(new SpectrumFile<TSimType, typename MagneticDiffusionTraits<TSimType>::MagType>(this->mMagB, "mag"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pMagSpectrum);
   }

}

#endif // MAGNETICDIFFUSIONSIMULATION_HPP
