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
#include "Simulations/SimulationBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a file generator
    *
    * \tparam TSimType Type of the simulation
    * \tparam TGenTraits Traits for the generator
    */
   template <typename TSimType, template <typename> class TGenTraits> class GeneratorBase: public SimulationBase<TSimType>
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
          * @brief Transform the field to spectral space
          */
         void transformRTP();

         /**
          * @brief Initialise the output file
          */
         virtual void initOutput(std::string name) = 0;

         /**
          * @brief write data to file
          */
         void writeFile();

         /**
          * @brief Finalise generator
          */
         void finalise();

         /**
          * @brief Get the codensity
          */
         typename TGenTraits<TSimType>::CodType&   codC();

         /**
          * @brief Get the codensity
          */
         typename TGenTraits<TSimType>::MagType&   magB();

         /**
          * @brief Get the codensity
          */
         typename TGenTraits<TSimType>::VelType&   velV();

      protected:
         /**
          * @brief Initialise the fields
          */
         void initFields();

         /**
          * @brief Initialise the output file
          *
          * @param pFile Output file smart pointer
          */
         void initOutputFile(SmartHDF5Writer pFile);

         /**
          * @brief Configure the transform manipulators
          */
         void configureTransform();

      private:

         /**
          * @brief Smart pointer of Codensity type
          */
         EPMSHARED_PTR<typename TGenTraits<TSimType>::CodType> mpCodC;

         /**
          * @brief Smart pointer of Codensity type
          */
         EPMSHARED_PTR<typename TGenTraits<TSimType>::MagType> mpMagB;

         /**
          * @brief Smart pointer of Codensity type
          */
         EPMSHARED_PTR<typename TGenTraits<TSimType>::VelType> mpVelV;
   };

   template <typename TSimType, template <typename> class TGenTraits> inline typename TGenTraits<TSimType>::CodType& GeneratorBase<TSimType, TGenTraits>::codC()
   {
      // Protect call with an assert
      assert(this->mpCodC != NULL);

      return (*this->mpCodC);
   }

   template <typename TSimType, template <typename> class TGenTraits> inline typename TGenTraits<TSimType>::MagType& GeneratorBase<TSimType, TGenTraits>::magB()
   {
      // Protect call with an assert
      assert(this->mpMagB != NULL);

      return (*this->mpMagB);
   }

   template <typename TSimType, template <typename> class TGenTraits> inline typename TGenTraits<TSimType>::VelType& GeneratorBase<TSimType, TGenTraits>::velV()
   {
      // Protect call with an assert
      assert(this->mpVelV != NULL);

      return (*this->mpVelV);
   }

   template <typename TSimType, template <typename> class TGenTraits> GeneratorBase<TSimType, TGenTraits>::GeneratorBase()
      : SimulationBase<TSimType>()
   {
      // Initialise fields
      this->initFields();

      // Configure transform
      this->configureTransform();
   }

   template <typename TSimType, template <typename> class TGenTraits> void GeneratorBase<TSimType, TGenTraits>::initFields()
   {
      // Create codensity field storage
      if(TGenTraits<TSimType>::NeedCodensity)
      {
         this->mpCodC.reset(new typename TGenTraits<TSimType>::CodType(this->mpTrunc, this->mTransform));
      }

      // Create magnetic field storage
      if(TGenTraits<TSimType>::NeedMagnetic)
      {
         this->mpMagB.reset(new typename TGenTraits<TSimType>::MagType(this->mpTrunc, this->mTransform));
      }

      // Create velocity field storage
      if(TGenTraits<TSimType>::NeedVelocity)
      {
         this->mpVelV.reset(new typename TGenTraits<TSimType>::VelType(this->mpTrunc, this->mTransform));
      }
   }

   template <typename TSimType, template <typename> class TGenTraits> void GeneratorBase<TSimType, TGenTraits>::configureTransform()
   {
      //
      // Setup the SSH transform data manipulator
      //

      // Register transform data packs of codensity scalar for SSH
      if(TGenTraits<TSimType>::NeedCodensity)
      {
         this->registerSSHPacks(1, 1);
      }

      // Register transform data packs of magnetic scalar for SSH
      if(TGenTraits<TSimType>::NeedMagnetic)
      {
         this->registerSSHPacks(3, 3);
      }

      // Register transform data packs of velocity scalar for SSH
      if(TGenTraits<TSimType>::NeedVelocity)
      {
         this->registerSSHPacks(3, 3);
      }

      // Configure the SSH manipulator
      this->configureSSHManipulator();

      //
      // Setup the SH transform data manipulator
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

      // Register transform data packs of codensity scalar for SH
      if(TGenTraits<TSimType>::NeedCodensity)
      {
         this->registerSHPacks(1, 1);
      }

      // Register transform data packs of magnetic scalar for SH
      if(TGenTraits<TSimType>::NeedMagnetic)
      {
         this->registerSHPacks(nSHBPacks, nSHFPacks);
      }

      // Register transform data packs of velocity scalar for SH
      if(TGenTraits<TSimType>::NeedVelocity)
      {
         this->registerSHPacks(nSHBPacks, nSHFPacks);
      }

      // Configure the SSH manipulator
      this->configureSHManipulator();

      //
      // Configure transform nesting
      //

      // Configure the transforms' nesting setup
      this->configureTransformNesting();
   }

   template <typename TSimType, template <typename> class TGenTraits> void GeneratorBase<TSimType, TGenTraits>::transformRTP()
   {
      for(int i=0; i < this->mTransformSteps; ++i)
      {
         // Control the groupe communications if applicable
         this->combineSpectralTransforms(i);

         // Transform the codensity field
         if(TGenTraits<TSimType>::NeedCodensity)
         {
            this->mTransform.transformRTP2Spec(this->codC().rOc().rPerturbation(), this->codC().oc().rtp());
         }

         // Transform the codensity field
         if(TGenTraits<TSimType>::NeedMagnetic)
         {
            this->mTransform.transformRTP2TorPol(this->magB().rOc().rPerturbation(), this->magB().oc().rtp());
         }

         // Transform the codensity field
         if(TGenTraits<TSimType>::NeedVelocity)
         {
            this->mTransform.transformRTP2TorPol(this->velV().rOc().rPerturbation(), this->velV().oc().rtp());
         }

      }
   }

   template <typename TSimType, template <typename> class TGenTraits> void GeneratorBase<TSimType, TGenTraits>::initOutputFile(SmartHDF5Writer pFile)
   {
      // Add the given state file to the output
      this->mIOSys.addHDF5Writer(pFile);

      // Initialise the writer
      this->mIOSys.initWriters();
   }

   template <typename TSimType, template <typename> class TGenTraits> void GeneratorBase<TSimType, TGenTraits>::writeFile()
   {
      // Write the created state to HDF5 file
      this->mIOSys.writeHDF5();
   }

   template <typename TSimType, template <typename> class TGenTraits> void GeneratorBase<TSimType, TGenTraits>::finalise()
   {
      // Print timestepping infos
      this->mSimControl.printInfo();

      // Close and finalise writers
      this->mIOSys.finaliseWriters();
   }

}

#endif // GENERATORBASE_HPP
