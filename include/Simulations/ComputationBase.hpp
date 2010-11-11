/** \file ComputationBase.hpp
 *  \brief Computation type independent parts
 */

#ifndef COMPUTATIONBASE_HPP
#define COMPUTATIONBASE_HPP

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
#include "General/ExecutionTimer.hpp"
#include "Domain/Truncation.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * @brief Implements the implementation independent part of a simulation
    *
    * \tparam TSimType General simulation type
    * \tparam TIOType IOSystem type
    */
   template <typename TSimType, typename TIOType> class ComputationBase
   {
      public:
         /// Typedef for the transform type
         typedef  typename SimulationTraits<TSimType>::TransformType  TransformType;

         /// Typedef for the EquationParameters type
         typedef typename SimulationTraits<TSimType>::EquationParametersType EquationParametersType;

         /**
          * @brief Simple empty destructor
          */
         virtual ~ComputationBase() {};

      protected:
         /**
          * @brief General execution timer
          */
         ExecutionTimer mExecTimer;

         /**
          * @brief Constructor
          */
         ComputationBase();

         /**
          * @brief IOSystem of the computation
          */
         TIOType    mIOSys;

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
         EquationParametersType   mEqParams;

         /**
          * @brief Number of steps required for a transform
          */
         int mTransformSteps;

         /**
          * @brief Do operations required just before starting the computation
          */
         void preRun();

         /**
          * @brief Do operations required just after finishing the computation
          */
         void postRun();

         /**
          * @brief Register maximum Forward and maximum packs for SSH transform
          *
          * This routine computes the actual maximum packs depending on the different call to this routine. This will most likely simply be called once for each considered equation.
          */
         void registerSSHPacks(const int maxFPacks, const int maxBPacks);

         /**
          * @brief Register maximum Forward and maximum packs for SH transform
          *
          * This routine computes the actual maximum packs depending on the different call to this routine. This will most likely simply be called once for each considered equation.
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
          * \epmTodo The multistep transform approach needs a review and maybe a refactoring
          */
         void combineRTPTransforms(const int entry);

         /**
          * @brief Combine RTP 2 Spectral transforms
          *
          * \epmTodo The multistep transform approach needs a review and maybe a refactoring
          */
         void combineSpectralTransforms(const int entry);

         /**
         * @name Methods that need to be overloaded in implementation
         */
         //@{

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

   template <typename TSimType, typename TIOType> ComputationBase<TSimType, TIOType>::ComputationBase()
      : mExecTimer(true), mIOSys(EquationParametersType::id), mpTrunc(TSimType::createTrunc(mIOSys.cfg()->aTrunc())), mTransform(mpTrunc), mEqParams(mIOSys.cfg()->aEquation()), mTransformSteps(0), mSSHFPacks(0), mSSHBPacks(0), mSHFPacks(0), mSHBPacks(0)
   {
      // Finish initialisation of the truncation object by setting the physical grid values
      this->mTransform.initRTPDomains(this->mpTrunc);
   }

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::combineRTPTransforms(const int entry)
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

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::combineSpectralTransforms(const int entry)
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

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::registerSSHPacks(const int maxFPacks, const int maxBPacks)
   {
      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         this->mSSHFPacks += maxFPacks; 
         this->mSSHBPacks += maxBPacks; 
      #else
         this->mSSHFPacks = std::max(this->mSSHFPacks, maxFPacks); 
         this->mSSHBPacks = std::max(this->mSSHBPacks, maxBPacks); 
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
   }

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::registerSHPacks(const int maxFPacks, const int maxBPacks)
   {
      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         this->mSHFPacks += maxFPacks; 
         this->mSHBPacks += maxBPacks; 
      #else
         this->mSHFPacks = std::max(this->mSHFPacks, maxFPacks); 
         this->mSHBPacks = std::max(this->mSHBPacks, maxBPacks); 
      #endif // EPMDYNAMO_SH_GROUPEDCOMM
   }

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::configureSSHManipulator()
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

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::configureSHManipulator()
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

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::configureTransformNesting()
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

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::preRun()
   {
      // Synchronize CPUs
      EPMDYNAMO_SYNCHRONIZE;

      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::postRun()
   {
      // Synchronize CPUs
      EPMDYNAMO_SYNCHRONIZE;

      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

   template <typename TSimType, typename TIOType> void ComputationBase<TSimType, TIOType>::configureTransforms()
   {
      BOOST_STATIC_ASSERT(sizeof(TSimType) == 0); 
   }

}

#endif // COMPUTATIONBASE_HPP
