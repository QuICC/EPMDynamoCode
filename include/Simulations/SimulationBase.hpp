/** \file SimulationBase.hpp
 *  \brief Implementation independent parts of a simulation
 */

#ifndef SIMULATIONBASE_HPP
#define SIMULATIONBASE_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <boost/static_assert.hpp>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Simulations/SimulationControl.hpp"
#include "IO/IOSystem.hpp"

namespace EPMDynamo {

   /**
    * @brief Implements the implementation independent part of a simulation
    *
    * \tparam TSimType General simulation type
    */
   template <typename TSimType> class SimulationBase
   {
      public:
         /// Typedef for the transform type
         typedef  typename SimulationTraits<TSimType>::TransformType  TransformType;

         /**
          * @brief Simple empty destructor
          */
         virtual ~SimulationBase() {};

      protected:
         /**
          * @brief Constructor
          */
         SimulationBase();

         /**
          * @brief IOSystem of the simulation
          */
         IOSystem    mIOSys;

         /**
          * @brief Shared pointer to the truncation
          */
         SmartTruncation  mpTrunc;

         /**
          * @brief Real space <-> spectral space transform object
          */
         TransformType  mTransform;

         /**
          * @brief Equation parameters
          */
         EquationParameters   mEqParams;

         /**
          * @brief SimulationBase control
          */
         SimulationControl<TSimType>   mSimControl;

         /**
          * @brief Number of steps required for a transform
          */
         int mTransformSteps;

         /**
          * @brief write data to files
          *
          * \epmBug This methods needs to be rename and cleaned out
          */
         void writeFiles();

         /**
          * @brief Initialise the IO output
          */
         void initOutput();

         /**
          * @brief Do operations required just before starting the time integration
          */
         void preRun();

         /**
          * @brief Do operations required just after finishing the time integration
          */
         void postRun();

         /**
          * @brief Register maximum Forward and maximum packs for SSH transform
          *
          * This routine computes the actual maximum packs depending on the different call to this routine. This will most likely simply be call once for each considered equation.
          */
         void registerSSHPacks(const int maxFPacks, const int maxBPacks);

         /**
          * @brief Register maximum Forward and maximum packs for SH transform
          *
          * This routine computes the actual maximum packs depending on the different call to this routine. This will most likely simply be call once for each considered equation.
          */
         void registerSHPacks(const int maxFPacks, const int maxBPacks);

         /**
          * @brief Configure the SpectralSH data manipulator
          */
         void configureSSHManipulator();

         /**
          * @brief Configure the SH data manipulator
          */
         void configureSHManipulator();

         /**
          * @brief Configure the nesting of the transforms
          */
         void configureTransformNesting();

         /**
          * @brief Combine Spectral 2 RTP transforms
          *
          * \bug This methods needs some more thoughts
          */
         void combineRTPTransforms(const int entry);

         /**
          * @brief Combine RTP 2 Spectral transforms
          *
          * \bug This methods needs some more thoughts
          */
         void combineSpectralTransforms(const int entry);

         /**
         * @name Methods that need to be overloaded in implementation
         */
         //@{

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
          * \param step Step in the transform for multistep trransforms
          */
         void updateEquationsRTP(const int step);

         /**
          * @brief Update equations RHS
          */
         void updateEquationsRHS();

         /**
          * @brief Transform equations RHS
          *
          * \param step Step in the transform for multistep trransforms
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
         //@}

      private:
         /**
          * @brief Storage for the registered number of SSH forward packs
          */
         int mSSHFPacks;

         /**
          * @brief Storage for the registered number of SSH backward packs
          */
         int mSSHBPacks;

         /**
          * @brief Storage for the registered number of SH forward packs
          */
         int mSHFPacks;

         /**
          * @brief Storage for the registered number of SH backward packs
          */
         int mSHBPacks;
   };

   template <typename TSimType> SimulationBase<TSimType>::SimulationBase()
      : mIOSys(), mpTrunc(TSimType::createTrunc(mIOSys.aTrunc())), mTransform(mpTrunc), mEqParams(mIOSys.aEquation()), mSimControl(mIOSys.aTStep(), mEqParams, mIOSys.aRunI(), mIOSys.aRun()), mTransformSteps(0), mSSHFPacks(0), mSSHBPacks(0), mSHFPacks(0), mSHBPacks(0)
   {
      // Finish initialisation of the truncation object by setting the physical grid values
      this->mTransform.initRTPDomains(this->mpTrunc);
   }

   template <typename TSimType> void SimulationBase<TSimType>::initOutput()
   {
      // Initialise all the create writers
      this->mIOSys.initWriters();
   }

   template <typename TSimType> void SimulationBase<TSimType>::combineRTPTransforms(const int entry)
   {
      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         this->mTransform.sshManipulator().initiateGroupedBSend(entry);
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM

      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         #ifndef EPMDYNAMO_RADIAL_GROUPEDCOMM
            this->mTransform.sshManipulator().initiateGroupedBSend(entry);
         #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
         this->mTransform.shManipulator().initiateGroupedBSend(entry);
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimType> void SimulationBase<TSimType>::combineSpectralTransforms(const int entry)
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         #ifndef EPMDYNAMO_RADIAL_GROUPEDCOMM
            this->mTransform.sshManipulator().initiateGroupedFSend(entry);
         #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
         this->mTransform.shManipulator().initiateGroupedFSend(entry);
      #endif // EPMDYNAMO_SHGROUPEDCOMM

      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         this->mTransform.sshManipulator().initiateGroupedFSend(entry);
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
   }

   template <typename TSimType> void SimulationBase<TSimType>::registerSSHPacks(const int maxFPacks, const int maxBPacks)
   {
      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         this->mSSHFPacks += maxFPacks; 
         this->mSSHBPacks += maxBPacks; 
      #else
         this->mSSHFPacks = std::max(this->mSSHFPacks, maxFPacks); 
         this->mSSHBPacks = std::max(this->mSSHBPacks, maxBPacks); 
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
   }

   template <typename TSimType> void SimulationBase<TSimType>::registerSHPacks(const int maxFPacks, const int maxBPacks)
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         this->mSHFPacks += maxFPacks; 
         this->mSHBPacks += maxBPacks; 
      #else
         this->mSHFPacks = std::max(this->mSHFPacks, maxFPacks); 
         this->mSHBPacks = std::max(this->mSHBPacks, maxBPacks); 
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimType> void SimulationBase<TSimType>::configureSSHManipulator()
   {
      bool state;

      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         state = false;
         this->mTransformSteps = 2;
      #else
         state = true;
         this->mTransformSteps = 1;
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM

      // Set the number of packs
      this->mTransform.sshManipulator().setMaxPacks(this->mSSHFPacks, this->mSSHBPacks);

      // Desactivate automatic send/recv control
      this->mTransform.sshManipulator().setDynamicPacks(state);

      // Activate or desactivate the special communcation entries
      this->mTransform.sshManipulator().setDesactivator(state);

      // Finish initialisation of manipulator
      this->mTransform.sshManipulator().setup();
   }

   template <typename TSimType> void SimulationBase<TSimType>::configureSHManipulator()
   {
      bool state;

      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         state = false;
         this->mTransformSteps++;
      #else
         state = true;
         this->mTransformSteps = std::max(this->mTransformSteps, 1);
      #endif // EPMDYNAMO_SH_GROUPEDCOMM

      // Set the number of packs
      this->mTransform.shManipulator().setMaxPacks(this->mSHFPacks, this->mSHBPacks);

      // Desactivate automatic send/recv control
      this->mTransform.shManipulator().setDynamicPacks(state);

      // Activate or desactivate the special communcation entries
      this->mTransform.shManipulator().setDesactivator(state);

      // Finish initialisation of manipulator
      this->mTransform.shManipulator().setup();
   }

   template <typename TSimType> void SimulationBase<TSimType>::configureTransformNesting()
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         // Set the SpectralSH manipulator to be aware of nested grouped communication
         this->mTransform.sshManipulator().setInterStageIDs(0, this->mTransformSteps-1);
         #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
            // If its a full massive grouped communication introduce index shift
            this->mTransform.shManipulator().setEntryShift(1);
         #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimType> void SimulationBase<TSimType>::writeFiles()
   {
      if(this->mSimControl.tsParams().isNextStep())
      {
         // Save ASCII files
         if(this->mSimControl.tsCounter().triggerASCII())
         {
            this->mIOSys.writeASCII();

            this->mSimControl.checkInterface();
         }

         // Save state file
         if(this->mSimControl.tsCounter().triggerState())
         {
            this->mIOSys.writeHDF5();
         }
      }
   }

   template <typename TSimType> void SimulationBase<TSimType>::preRun()
   {
      // Write loaded state to HDF5 file (this can be important if the truncations were not the same)
      this->mIOSys.writeHDF5();

      // Synchronize simulation
      EPMDYNAMO_SYNCHRONIZE;

      // Execute last minute simulation control initialisation operations
      this->mSimControl.preRun();
   }

   template <typename TSimType> void SimulationBase<TSimType>::postRun()
   {
      // Write the state to HDF5 file before finishing
      this->mIOSys.writeHDF5();

      // Execute early simulation control finalisation operations
      this->mSimControl.postRun();
   }

   template <typename TSimType> void SimulationBase<TSimType>::initFields()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::addHDF5Output()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::addASCIIOutput()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::initEquations()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::updateEquationsRTP(const int step)
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::updateEquationsRHS()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::transformEquationsRHS(const int step)
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::addExternalInfluence()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::timestepEquations()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType> void SimulationBase<TSimType>::configureTransforms()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

}

#endif // SIMULATIONBASE_HPP
