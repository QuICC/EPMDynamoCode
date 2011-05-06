/** \file InitialStateGenerator.hpp
 *  \brief Implementation of a initial state generator
 */

#ifndef INITIALSTATEGENERATOR_HPP
#define INITIALSTATEGENERATOR_HPP

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
#include "Utilities/GeneratorBase.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a few methods for easier implementation of an initial state generator
    *
    * \tparam TSimTraits Traits for the generator
    */
   template <typename TSimTraits> class InitialStateGenerator: public GeneratorBase<TSimTraits>
   {
      public:
         /// Typedef for a smart state file
         typedef EPMSHARED_PTR<StateFileWriter<TSimTraits> > SmartStateWriter; 

         /**
          * @brief Constructor
          */
         InitialStateGenerator();

         /**
          * @brief Simple empty destructor
          */
         virtual ~InitialStateGenerator() {};

         /**
          * @brief Initialise the output file
          *
          * @param name Basename of the output file
          */
         void setupOutput(std::string name = "Initial");

         /**
          * @brief Write the state file
          */
         void writeStateFile();

         /**
          * @brief Finalise the output state file
          */
         void finalise();

      protected:
         /**
          * @brief Smart pointer for HDF5 output
          */
         SmartStateWriter   mpOutFile;

         /**
          * @brief Initialise the output state file
          *
          * @param pFile Output file smart pointer
          */
         void initStateFile(SmartStateWriter pFile);


      private:
   };

   template <typename TSimTraits> InitialStateGenerator<TSimTraits>::InitialStateGenerator()
      : GeneratorBase<TSimTraits>()
   {
   }

   template <typename TSimTraits> void InitialStateGenerator<TSimTraits>::setupOutput(std::string name)
   {
      SmartStateWriter  pOutFile;

      // Create state with all fields
      if(TSimTraits::NeedCodensity && TSimTraits::NeedMagnetic && TSimTraits::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimTraits>(this->codC(), this->magB(), this->velV(), this->mEqParams, this->mTSParams));

      // Create state with Codensity and Velocity fields
      } else if(TSimTraits::NeedCodensity && TSimTraits::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimTraits>(this->codC(), this->velV(), this->mEqParams, this->mTSParams));

      // Create state with Magnetic and Velocity fields
      } else if(TSimTraits::NeedMagnetic && TSimTraits::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimTraits>(this->magB(), this->velV(), this->mEqParams, this->mTSParams));

      // Create state with only codensity field
      } else if(TSimTraits::NeedCodensity)
      {
         pOutFile.reset(new StateFileWriter<TSimTraits>(this->codC(), this->mEqParams, this->mTSParams));

      // Create state with only magnetic field
      } else if(TSimTraits::NeedMagnetic)
      {
         pOutFile.reset(new StateFileWriter<TSimTraits>(this->magB(), this->mEqParams, this->mTSParams));

      // Create state with only velocity field
      } else if(TSimTraits::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimTraits>(this->velV(), this->mEqParams, this->mTSParams));
      }

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initStateFile(pOutFile);
   }

   template <typename TSimTraits> void InitialStateGenerator<TSimTraits>::initStateFile(typename InitialStateGenerator<TSimTraits>::SmartStateWriter pFile)
   {
      // Set output file
      this->mpOutFile = pFile;
   }

   template <typename TSimTraits> void InitialStateGenerator<TSimTraits>::writeStateFile()
   {
      this->mpOutFile->write();
   }

   template <typename TSimTraits> void InitialStateGenerator<TSimTraits>::finalise()
   {
      // Finalise output file
      this->mpOutFile->finalise();

      // Finalise generator
      GeneratorBase<TSimTraits>::finalise();
   }

}

#endif // INITIALSTATEGENERATOR_HPP
