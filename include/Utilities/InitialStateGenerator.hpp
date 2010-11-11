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
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits for the generator
    */
   template <typename TSimType, template <typename> class TSimTraits> class InitialStateGenerator: public GeneratorBase<TSimType, TSimTraits>
   {
      public:
         /// Typedef for a smart state file
         typedef EPMSHARED_PTR<StateFileWriter<TSimType, TSimTraits> > SmartStateWriter; 

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

   template <typename TSimType, template <typename> class TSimTraits> InitialStateGenerator<TSimType, TSimTraits>::InitialStateGenerator()
      : GeneratorBase<TSimType, TSimTraits>()
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void InitialStateGenerator<TSimType, TSimTraits>::setupOutput(std::string name)
   {
      SmartStateWriter  pOutFile;

      // Create state with all fields
      if(TSimTraits<TSimType>::NeedCodensity && TSimTraits<TSimType>::NeedMagnetic && TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TSimTraits>(this->codC(), this->magB(), this->velV(), this->mEqParams, this->mTSParams));

      // Create state with Codensity and Velocity fields
      } else if(TSimTraits<TSimType>::NeedCodensity && TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TSimTraits>(this->codC(), this->velV(), this->mEqParams, this->mTSParams));

      // Create state with Magnetic and Velocity fields
      } else if(TSimTraits<TSimType>::NeedMagnetic && TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TSimTraits>(this->magB(), this->velV(), this->mEqParams, this->mTSParams));

      // Create state with only codensity field
      } else if(TSimTraits<TSimType>::NeedCodensity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TSimTraits>(this->codC(), this->mEqParams, this->mTSParams));

      // Create state with only magnetic field
      } else if(TSimTraits<TSimType>::NeedMagnetic)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TSimTraits>(this->magB(), this->mEqParams, this->mTSParams));

      // Create state with only velocity field
      } else if(TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TSimTraits>(this->velV(), this->mEqParams, this->mTSParams));
      }

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initStateFile(pOutFile);
   }

   template <typename TSimType, template <typename> class TSimTraits> void InitialStateGenerator<TSimType, TSimTraits>::initStateFile(typename InitialStateGenerator<TSimType, TSimTraits>::SmartStateWriter pFile)
   {
      // Set output file
      this->mpOutFile = pFile;
   }

   template <typename TSimType, template <typename> class TSimTraits> void InitialStateGenerator<TSimType, TSimTraits>::writeStateFile()
   {
      this->mpOutFile->write();
   }

   template <typename TSimType, template <typename> class TSimTraits> void InitialStateGenerator<TSimType, TSimTraits>::finalise()
   {
      // Finalise output file
      this->mpOutFile->finalise();

      // Finalise generator
      GeneratorBase<TSimType, TSimTraits>::finalise();
   }

}

#endif // INITIALSTATEGENERATOR_HPP
