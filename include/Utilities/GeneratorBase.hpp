/** \file GeneratorBase.hpp
 *  \brief Base of the implementation of a file generator
 */

#ifndef GENERATORBASE_HPP
#define GENERATORBASE_HPP

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
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Simulations/ComputationBase.hpp"
#include "IO/IOBaseSystem.hpp"
#include "IO/HDF5/HDF5Reader.hpp"
#include "IO/HDF5/HDF5Writer.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a file generator
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits for the generator
    */
   template <typename TSimType, template <typename> class TSimTraits> class GeneratorBase: public ComputationBase<TSimType, IOBaseSystem>
   {
      public:
         /// Typedef for the transform type
         typedef  typename SimulationTraits<TSimType>::TransformType   TransformType;

         /**
          * @brief Constructor
          */
         GeneratorBase();

         /**
          * @brief Simple empty destructor
          */
         virtual ~GeneratorBase() {};

         /**
          * @brief Get the codensity
          */
         typename TSimTraits<TSimType>::CodType&   codC();

         /**
          * @brief Get the codensity
          */
         typename TSimTraits<TSimType>::MagType&   magB();

         /**
          * @brief Get the codensity
          */
         typename TSimTraits<TSimType>::VelType&   velV();

         /**
          * @brief Configure the transform manipulators
          */
         template <typename TGenTraits> void configureTransform();

         /**
          * @brief Transform the (required) fields to spectral space
          */
         template <typename TGenTraits> void transformRTP();

         /**
          * @brief Transform the (required) fields to RTP space
          */
         template <typename TGenTraits> void transformSpectral();

         /**
          * @brief Finalise generator
          */
         void finalise();

      protected:
         /**
          * @brief Initialise the (required) fields
          */
         void initFields();

         /**
          * @brief Timestep parameters
          */
         TimestepParameters   mTSParams;

      private:
         /**
          * @brief Smart pointer of Codensity type
          */
         EPMSHARED_PTR<typename TSimTraits<TSimType>::CodType> mpCodC;

         /**
          * @brief Smart pointer of Codensity type
          */
         EPMSHARED_PTR<typename TSimTraits<TSimType>::MagType> mpMagB;

         /**
          * @brief Smart pointer of Codensity type
          */
         EPMSHARED_PTR<typename TSimTraits<TSimType>::VelType> mpVelV;


   };

   template <typename TSimType, template <typename> class TSimTraits> inline typename TSimTraits<TSimType>::CodType& GeneratorBase<TSimType, TSimTraits>::codC()
   {
      // Protect call with an assert
      assert(this->mpCodC != NULL);

      return (*this->mpCodC);
   }

   template <typename TSimType, template <typename> class TSimTraits> inline typename TSimTraits<TSimType>::MagType& GeneratorBase<TSimType, TSimTraits>::magB()
   {
      // Protect call with an assert
      assert(this->mpMagB != NULL);

      return (*this->mpMagB);
   }

   template <typename TSimType, template <typename> class TSimTraits> inline typename TSimTraits<TSimType>::VelType& GeneratorBase<TSimType, TSimTraits>::velV()
   {
      // Protect call with an assert
      assert(this->mpVelV != NULL);

      return (*this->mpVelV);
   }



   template <typename TSimType, template <typename> class TSimTraits> GeneratorBase<TSimType, TSimTraits>::GeneratorBase()
      : ComputationBase<TSimType, IOBaseSystem>(), mTSParams(this->mIOSys.cfg()->aTStep()(0), this->mIOSys.cfg()->aTStep()(1))
   {
      // Initialise fields
      this->initFields();
   }

   template <typename TSimType, template <typename> class TSimTraits> void GeneratorBase<TSimType, TSimTraits>::initFields()
   {
      // Create codensity field storage
      if(TSimTraits<TSimType>::NeedCodensity)
      {
         this->mpCodC.reset(new typename TSimTraits<TSimType>::CodType(this->mpTrunc, this->mTransform));
      }

      // Create magnetic field storage
      if(TSimTraits<TSimType>::NeedMagnetic)
      {
         this->mpMagB.reset(new typename TSimTraits<TSimType>::MagType(this->mpTrunc, this->mTransform));
      }

      // Create velocity field storage
      if(TSimTraits<TSimType>::NeedVelocity)
      {
         this->mpVelV.reset(new typename TSimTraits<TSimType>::VelType(this->mpTrunc, this->mTransform));
      }
   }

   template <typename TSimType, template <typename> class TSimTraits> template <typename TGenTraits> void GeneratorBase<TSimType, TSimTraits>::configureTransform()
   {
      //
      // Setup the SSH and SH transform data manipulator
      //

      int nSHFPacks;
      int nSHBPacks;

      #ifdef EPMDYNAMO_SH_GROUPEDCOMM
         nSHFPacks = 3;
         nSHBPacks = 3;
      #else
         nSHFPacks = 2;
         nSHBPacks = 1;
      #endif // EPMDYNAMO_SH_GROUPEDCOMM

      // Register transform data packs of codensity scalar for SSH
      if(TSimTraits<TSimType>::NeedCodensity && TGenTraits::UseRTPCodensity)
      {
         // Register SSH packs
         this->registerSSHPacks(1, 0);

         // Register SH packs
         this->registerSHPacks(1, 0);
      }

      // Register transform data packs of codensity scalar for SSH
      if(TSimTraits<TSimType>::NeedCodensity && TGenTraits::UseSpecCodensity)
      {
         // Register SSH packs
         this->registerSSHPacks(0, 1);

         // Register SH packs
         this->registerSHPacks(0, 1);
      }

      // Register transform data packs of codensity gradient for SSH
      if(TSimTraits<TSimType>::NeedCodensity && TGenTraits::UseCodensityGrad)
      {
         // Register SSH packs
         this->registerSSHPacks(0, 1);

         // Register SH packs
         this->registerSHPacks(0, 1);
      }

      // Register transform data packs of magnetic field for SSH
      if(TSimTraits<TSimType>::NeedMagnetic && TGenTraits::UseRTPMagnetic)
      {
         // Register SSH packs
         this->registerSSHPacks(3, 0);

         // Register SH packs
         this->registerSHPacks(nSHFPacks, 0);
      }

      // Register transform data packs of magnetic field for SSH
      if(TSimTraits<TSimType>::NeedMagnetic && TGenTraits::UseSpecMagnetic)
      {
         // Register SSH packs
         this->registerSSHPacks(0, 3);

         // Register SH packs
         this->registerSHPacks(0, nSHBPacks);
      }

      // Register transform data packs of magnetic curl for SSH
      if(TSimTraits<TSimType>::NeedMagnetic && TGenTraits::UseMagneticCurl)
      {
         // Register SSH packs
         this->registerSSHPacks(0, 3);

         // Register SH packs
         this->registerSHPacks(0, nSHBPacks);
      }

      // Register transform data packs of velocity field for SSH
      if(TSimTraits<TSimType>::NeedVelocity && TGenTraits::UseRTPVelocity)
      {
         // Register SSH packs
         this->registerSSHPacks(3, 0);

         // Register SH packs
         this->registerSHPacks(nSHFPacks, 0);
      }

      // Register transform data packs of velocity field for SSH
      if(TSimTraits<TSimType>::NeedVelocity && TGenTraits::UseSpecVelocity)
      {
         // Register SSH packs
         this->registerSSHPacks(0, 3);

         // Register SH packs
         this->registerSHPacks(0, nSHBPacks);
      }

      // Register transform data packs of velocity curl for SSH
      if(TSimTraits<TSimType>::NeedVelocity && TGenTraits::UseVelocityCurl)
      {
         // Register SSH packs
         this->registerSSHPacks(0, 3);

         // Register SH packs
         this->registerSHPacks(0, nSHBPacks);
      }

      // Configure the SSH manipulator
      this->configureSSHManipulator();

      // Configure the SH manipulator
      this->configureSHManipulator();

      //
      // Configure transform nesting
      //

      // Configure the transforms' nesting setup
      this->configureTransformNesting();
   }

   template <typename TSimType, template <typename> class TSimTraits>  template <typename TGenTraits> void GeneratorBase<TSimType, TSimTraits>::transformRTP()
   {
      for(int i=0; i < this->mTransformSteps; ++i)
      {
         // Control the group communications if applicable
         this->combineSpectralTransforms(i);

         // Transform the codensity field
         if(TSimTraits<TSimType>::NeedCodensity && TGenTraits::UseRTPCodensity)
         {
            this->mTransform.transformRTP2Spec(this->codC().rOc().rPerturbation(), this->codC().oc().rtp());
         }

         // Transform the codensity field
         if(TSimTraits<TSimType>::NeedMagnetic && TGenTraits::UseRTPMagnetic)
         {
            this->mTransform.transformRTP2TorPol(this->magB().rOc().rPerturbation(), this->magB().oc().rtp());
         }

         // Transform the codensity field
         if(TSimTraits<TSimType>::NeedVelocity && TGenTraits::UseRTPVelocity)
         {
            this->mTransform.transformRTP2TorPol(this->velV().rOc().rPerturbation(), this->velV().oc().rtp());
         }

      }
   }

   template <typename TSimType, template <typename> class TSimTraits>  template <typename TGenTraits> void GeneratorBase<TSimType, TSimTraits>::transformSpectral()
   {
      for(int i=0; i < this->mTransformSteps; ++i)
      {
         // Control the groupe communications if applicable
         this->combineRTPTransforms(i);

         // Transform the codensity scalar
         if(TSimTraits<TSimType>::NeedCodensity && TGenTraits::UseSpecCodensity)
         {
            this->codC().rOc().transform(i);
         }

         // Transform the codensity gradient
         if(TSimTraits<TSimType>::NeedCodensity && TGenTraits::UseCodensityGrad)
         {
            this->codC().rOc().gradTransform(i);
         }

         // Transform the magnetic field
         if(TSimTraits<TSimType>::NeedMagnetic && TGenTraits::UseSpecMagnetic)
         {
            this->magB().rOc().transform(i);
         }

         // Transform the magnetic curl
         if(TSimTraits<TSimType>::NeedMagnetic && TGenTraits::UseMagneticCurl)
         {
            this->magB().rOc().curlTransform(i);
         }

         // Transform the velocity field
         if(TSimTraits<TSimType>::NeedVelocity && TGenTraits::UseSpecVelocity)
         {
            this->velV().rOc().transform(i);
         }

         // Transform the velocity curl
         if(TSimTraits<TSimType>::NeedVelocity && TGenTraits::UseVelocityCurl)
         {
            this->velV().rOc().curlTransform(i);
         }
      }
   }

   template <typename TSimType, template <typename> class TSimTraits> void GeneratorBase<TSimType, TSimTraits>::finalise()
   {
      // Finalise IOBaseSystem
      this->mIOSys.finaliseBase();
   }

}

#endif // GENERATORBASE_HPP
