/** \file SimulationBase.hpp
 *  \brief Implementation independent parts of a simulation
 */

#ifndef SIMULATIONBASE_HPP
#define SIMULATIONBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/SimulationTruncation.hpp"
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

         /// Typedef for the equation parameters type
         typedef  typename TSimType::EqParamsType   EqParamsType;

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
          * @brief Shared pointer to the simulation truncation
          */
         SmartSTrunc  mpSTrunc;

         /**
          * @brief Real space <-> spectral space transform object
          */
         TransformType  mTransform;

         /**
          * @brief Equation parameters
          */
         EqParamsType   mEqParams;

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
          * \bug This methods needs to be rename and cleaned out
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
          * @brief Configure the WSH data manipulator
          *
          * \bug This methods needs to be rename and cleaned out
          */
         void configureWSHManipulator(const int maxFPacks, const int maxBPacks, bool state);

         /**
          * @brief Configure the SH data manipulator
          *
          * \bug This methods needs to be rename and cleaned out
          */
         void configureSHManipulator(const int maxFPacks, const int maxBPacks, bool state);

         /**
          * @brief Configure the nested data manipulator
          *
          * \bug This methods needs to be rename and cleaned out
          */
         void configureNestedManipulators();

         /**
          * @brief Combine Spectral 2 RTP transforms
          *
          * \bug This methods needs to be rename and cleaned out
          */
         void combineRTPTransforms(const int entry);

         /**
          * @brief Combine RTP 2 Spectral transforms
          *
          * \bug This methods needs to be rename and cleaned out
          */
         void combineSpectralTransforms(const int entry);

      private:
   };

   template <typename TSimType> SimulationBase<TSimType>::SimulationBase()
      : mIOSys(), mpSTrunc(TSim::createSTrunc(mIOSys.aTrunc())), mTransform(mpSTrunc), mEqParams(mIOSys.aEquation()), mSimControl(mIOSys.aTStep(), mEqParams, mIOSys.aRunI(), mIOSys.aRun()), mTransformSteps(0)
   {
      // Finish initialisation of the truncation object by setting the physical grid values
      this->mTransform.initRTPDomains(this->mpSTrunc);
   }

   template <typename TSimType> void SimulationBase<TSimType>::initOutput()
   {
      // Initialise all the create writers
      this->mIOSys.initWriters();
   }

   template <typename TSimType> void SimulationBase<TSimType>::combineRTPTransforms(const int entry)
   {
      #ifdef EPMDYNAMO_GROUPED_ONESTEP
         this->mTransform.spectralSHManipulator().initiateGroupedBSend(entry);
      #endif // EPMDYNAMO_GROUPED_ONESTEP

      #ifdef EPMDYNAMO_GROUPED_TWOSTEP
         #ifndef EPMDYNAMO_GROUPED_ONESTEP
            this->mTransform.spectralSHManipulator().initiateGroupedBSend(entry);
         #endif // EPMDYNAMO_GROUPED_ONESTEP
         this->mTransform.shManipulator().initiateGroupedBSend(entry);
      #endif // EPMDYNAMO_GROUPED_TWOSTEP
   }

   template <typename TSimType> void SimulationBase<TSimType>::combineSpectralTransforms(const int entry)
   {
      #ifdef EPMDYNAMO_GROUPED_TWOSTEP
         #ifndef EPMDYNAMO_GROUPED_ONESTEP
            this->mTransform.spectralSHManipulator().initiateGroupedFSend(entry);
         #endif // EPMDYNAMO_GROUPED_ONESTEP
         this->mTransform.shManipulator().initiateGroupedFSend(entry);
      #endif // EPMDYNAMO_GROUPED_TWOSTEP

      #ifdef EPMDYNAMO_GROUPED_ONESTEP
         this->mTransform.spectralSHManipulator().initiateGroupedFSend(entry);
      #endif // EPMDYNAMO_GROUPED_ONESTEP
   }

   template <typename TSimType> void SimulationBase<TSimType>::configureWSHManipulator(const int maxFPacks, const int maxBPacks, bool state)
   {
      // Set the number of packs
      mTransform.spectralSHManipulator().setMaxPacks(maxFPacks, maxBPacks);

      // Desactivate automatic send/recv control
      mTransform.spectralSHManipulator().setDynamicPacks(state);

      // Activate or desactivate the special communcation entries
      mTransform.spectralSHManipulator().setDesactivator(state);

      // Finish initialisation of manipulator
      mTransform.spectralSHManipulator().setup();
   }

   template <typename TSimType> void SimulationBase<TSimType>::configureSHManipulator(const int maxFPacks, const int maxBPacks, bool state)
   {
      #ifdef EPMDYNAMO_TWOSTEP
         // Set the number of packs
         mTransform.shManipulator().setMaxPacks(maxFPacks, maxBPacks);

         // Desactivate automatic send/recv control
         mTransform.shManipulator().setDynamicPacks(state);

         // Activate or desactivate the special communcation entries
         mTransform.shManipulator().setDesactivator(state);

         // Finish initialisation of manipulator
         mTransform.shManipulator().setup();
      #endif // EPMDYNAMO_TWOSTEP
   }

   template <typename TSimType> void SimulationBase<TSimType>::configureNestedManipulators()
   {
      #ifdef EPMDYNAMO_GROUPED_TWOSTEP
         // Set the WSH manipulator to be aware of nested grouped communication
         mTransform.spectralSHManipulator().setSpecialEntries(0, this->mTransformSteps-1);
         #ifdef EPMDYNAMO_GROUPED_ONESTEP
            // If its a full massive grouped communication introduce index shift
            mTransform.shManipulator().setEntryShift(1);
         #endif // EPMDYNAMO_GROUPED_ONESTEP
      #endif // EPMDYNAMO_GROUPED_TWOSTEP
   }

   template <typename TSimType> void SimulationBase<TSimType>::writeFiles()
   {
      if(this->mSimControl.tsParams().isNextStep())
      {
         // Save ASCII files
         if(this->mSimControl.tsCounter().triggerA())
         {
            this->mIOSys.writeASCII();

            this->mSimControl.checkInterface();
         }

         // Save state file
         if(this->mSimControl.tsCounter().triggerS())
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
      this->mSimControl.synchronize();

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

}

#endif // SIMULATIONBASE_HPP
