/** \file CodensityDiffusionSimulation.hpp
 *  \brief Implementation of a codensity diffusion simulation
 */

#ifndef CODENSITYDIFFUSIONSIMULATION_HPP
#define CODENSITYDIFFUSIONSIMULATION_HPP

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
#include "Simulations/Traits/CodensityDiffusionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

#include "Equations/Transport/TransportDiffusion.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a codensity diffusion simulation
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class CodensityDiffusionSimulation: public SimulationBase<TSimType>
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~CodensityDiffusionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         CodensityDiffusionSimulation();

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
         typename CodensityDiffusionTraits<TSimType>::CodType   mCodC;

         /**
          * @brief Transport equation
          */
         TransportDiffusion<TSimType, CodensityDiffusionTraits>    mTransport;
   };

   template <typename TSimType> CodensityDiffusionSimulation<TSimType>::CodensityDiffusionSimulation()
      : mCFLFactor(1.0), mCodC(this->mpTrunc, this->mTransform), mTransport(mCodC, this->mTransform, this->mSimControl.tsParams(), this->mEqParams)
   {
      // Set the CFL factor to L*(L+1)
      int l = this->mpTrunc->sim()->hoz()->nL();
      this->mCFLFactor = static_cast<EPMFloat>(l*(l+1));
   }

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::initEquations()
   {
      // Set boundary condition to the transport equation
      SmartBC  pZeroBC(new ZeroBC<TSimType>(this->mTransform.radBasis()));
      this->mTransport.addBC(pZeroBC);

      // Initialise the transport equation
      this->mTransport.init();
   }

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::configureTransforms()
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

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::updateEquationsRTP(const int step)
   {
      // Update RTP values of the transport equation
      this->mTransport.updateRTP(step);
   }

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::updateEquationsRHS()
   {
      // Update RHS of the transport equation
      this->mTransport.updateRHS();
   }

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::transformEquationsRHS(const int step)
   {
      // Update RHS of the transport equation
      this->mTransport.transformRHS(step);
   }

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::addExternalInfluence()
   {
   }

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::timestepEquations()
   {
      // Timestep the transport equation
      this->mTransport.timestep();
   }

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::initFields()
   {
      // Create a state file reader for the initial state
      EPMSHARED_PTR<StateFileReader<TSimType, CodensityDiffusionTraits> > pInState(new StateFileReader<TSimType, CodensityDiffusionTraits>(this->mCodC,  "_initial"));

      // Read in initial state
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create a source file reader for the codensity source
      EPMSHARED_PTR<SourceFileReader<TSimType, CodensityDiffusionTraits> > pSource(new SourceFileReader<TSimType, CodensityDiffusionTraits>(this->mCodC));

      // Read in source state
      this->mIOSys.useSource(pSource);

      // Rescale the source term to correct parameter related value
      this->mCodC.rOcSrc().rescale(this->mEqParams.codSourceScale());
   }

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::addHDF5Output()
   {
      EPMSHARED_PTR<StateFileWriter<TSimType, CodensityDiffusionTraits> >  pOutState(new StateFileWriter<TSimType, CodensityDiffusionTraits>(this->mCodC, this->mEqParams, this->mSimControl.tsParams()));

      this->mIOSys.addHDF5Writer(pOutState);
   }

   template <typename TSimType> void CodensityDiffusionSimulation<TSimType>::addASCIIOutput()
   {
      // Create a timestep ASCII logging file
      EPMSHARED_PTR<TimeFile> pTimeFile(new TimeFile("timestep", this->mSimControl.tsParams()));
      // Add time file to ASCII output
      this->mIOSys.addASCIIWriter(pTimeFile);

      // Create a energy ASCII diagnostic file for the codensity field
      EPMSHARED_PTR<EnergyFile<TSimType, typename CodensityDiffusionTraits<TSimType>::CodType> > pCodEnergy(new EnergyFile<TSimType, typename CodensityDiffusionTraits<TSimType>::CodType>(mCodC, "cod", this->mSimControl.tsParams()));

      // Add kinetic energy to ASCII output
      this->mIOSys.addASCIIWriter(pCodEnergy);

      // Create a energy spectrum ASCII diagnostic file for the codensity field
      EPMSHARED_PTR<SpectrumFile<TSimType, typename CodensityDiffusionTraits<TSimType>::CodType> > pCodSpectrum(new SpectrumFile<TSimType, typename CodensityDiffusionTraits<TSimType>::CodType>(this->mCodC, "cod"));

      // Add kinetic energy spectrum to ASCII output
      this->mIOSys.addASCIIWriter(pCodSpectrum);
   }

}

#endif // CODENSITYDIFFUSIONSIMULATION_HPP
