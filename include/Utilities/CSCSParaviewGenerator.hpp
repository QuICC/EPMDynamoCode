/** \file CSCSParaviewGenerator.hpp
 *  \brief Implementation of a simple CSCS visualisation file format generator
 */

#ifndef CSCSPARAVIEWGENERATOR_HPP
#define CSCSPARAVIEWGENERATOR_HPP

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
#include "Utilities/Traits/VisGeneratorTraits.hpp"
#include "Utilities/Traits/VisStateFilterTraits.hpp"
#include "IO/HDF5/Visualisation/CSCSFileWriter.hpp"
#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/Visualisation/CSCSFileDefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a simple CSCS visualisation file format generator
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits for the generator
    */
   template <typename TSimType, template <typename> class TSimTraits> class CSCSParaviewGenerator: public GeneratorBase<TSimType, TSimTraits>
   {
      public:
         /// Typedef for the state file reader
         typedef EPMSHARED_PTR<StateFileReader<TSimType, TSimTraits> >  SmartStateReader; 

         /// Typedef for the CSCS file writer
         typedef EPMSHARED_PTR<CSCSFileWriter<TSimType, TSimTraits> >   SmartCSCSWriter;

         /**
          * @brief Constructor
          */
         CSCSParaviewGenerator();

         /**
          * @brief Simple empty destructor
          */
         virtual ~CSCSParaviewGenerator() {};

         /**
          * @brief Initialise the output file
          *
          * @param name Basename of the file
          */
         void setupOutput(std::string name = "CSCSParaView");

         /**
          * @brief Write total field values
          */
         template <typename TVisTraits> void writeTotal();

         /**
          * @brief Write toroidal component values
          */
         template <typename TVisTraits> void writeToroidal();

         /**
          * @brief Write poloidal component values
          */
         template <typename TVisTraits> void writePoloidal();

         /**
          * @brief Finalise the output cscs file
          */
         void finalise();

      protected:

         /**
          * @brief Smart pointer for the state file reader
          */
         SmartStateReader  mpInFile;

         /**
          * @brief Smart pointer for the CSCS file writer
          */
         SmartCSCSWriter  mpOutFile;

         /**
          * @brief Initialise the input file
          *
          * @param pFile Smart pointer to file
          */
         void initInputFile(SmartStateReader  pFile);

         /**
          * @brief Initialise the cscs output file
          *
          * @param pFile Smart pointer to file
          */
         void initCSCSFile(SmartCSCSWriter  pFile);

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> CSCSParaviewGenerator<TSimType, TSimTraits>::CSCSParaviewGenerator()
      : GeneratorBase<TSimType, TSimTraits>()
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void CSCSParaviewGenerator<TSimType, TSimTraits>::setupOutput(std::string name)
   {
      SmartStateReader  pInFile;
      SmartCSCSWriter  pOutFile;

      // Create file with all fields
      if(TSimTraits<TSimType>::NeedCodensity && TSimTraits<TSimType>::NeedMagnetic && TSimTraits<TSimType>::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimType, TSimTraits>(this->codC(), this->magB(), this->velV(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimType, TSimTraits>("FSOuter", this->codC(), this->magB(), this->velV(), this->mTSParams));

      // Create file with Codensity and Velocity fields
      } else if(TSimTraits<TSimType>::NeedCodensity && TSimTraits<TSimType>::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimType, TSimTraits>(this->codC(), this->velV(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimType, TSimTraits>("FSOuter", this->codC(), this->velV(), this->mTSParams));

      // Create file with Magnetic and Velocity fields
      } else if(TSimTraits<TSimType>::NeedMagnetic && TSimTraits<TSimType>::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimType, TSimTraits>(this->magB(), this->velV(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimType, TSimTraits>("FSOuter", this->magB(), this->velV(), this->mTSParams));

      // Create file with only codensity field
      } else if(TSimTraits<TSimType>::NeedCodensity)
      {
         pInFile.reset(new StateFileReader<TSimType, TSimTraits>(this->codC(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimType, TSimTraits>("FSOuter", this->codC(), this->mTSParams));

      // Create file with only magnetic field
      } else if(TSimTraits<TSimType>::NeedMagnetic)
      {
         pInFile.reset(new StateFileReader<TSimType, TSimTraits>(this->magB(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimType, TSimTraits>("FSOuter", this->magB(), this->mTSParams));

      // Create file with only velocity field
      } else if(TSimTraits<TSimType>::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimType, TSimTraits>(this->velV(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimType, TSimTraits>("FSOuter", this->velV(), this->mTSParams));
      }

      // Init the input file
      this->initInputFile(pInFile);

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initCSCSFile(pOutFile);
   }

   template <typename TSimType, template <typename> class TSimTraits> void CSCSParaviewGenerator<TSimType, TSimTraits>::initInputFile(typename CSCSParaviewGenerator<TSimType, TSimTraits>::SmartStateReader   pInFile)
   {
      // Set the state file localy
      this->mpInFile = pInFile;

      // Initialise the input file
      this->mpInFile->init();

      // Read setup information from file
      this->mpInFile->readSetup();
   }

   template <typename TSimType, template <typename> class TSimTraits> void CSCSParaviewGenerator<TSimType, TSimTraits>::initCSCSFile(typename CSCSParaviewGenerator<TSimType, TSimTraits>::SmartCSCSWriter   pFile)
   {
      // Set output file
      this->mpOutFile = pFile;

      // Set output file
      this->mpOutFile->init();

      // Write start information to file
      this->mpOutFile->writeStart();
   }

   template <typename TSimType, template <typename> class TSimTraits> template <typename TVisTraits> void CSCSParaviewGenerator<TSimType, TSimTraits>::writeTotal()
   {
      // Read data from file
      this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::FullField> >();

      // Configure transforms
      this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

      // Transform the fields
      this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

      // Write data to file
      this->mpOutFile->template writeVisualisation<TVisTraits>();
   }

   template <typename TSimType, template <typename> class TSimTraits> template <typename TVisTraits> void CSCSParaviewGenerator<TSimType, TSimTraits>::writeToroidal()
   {
      // Read data from file
      this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::ToroidalOnly> >();

      // Configure transforms
      this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

      // Transform the fields
      this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

      // Write data to file
      this->mpOutFile->template writeVisualisation<TVisTraits>(CSCSFileDefs::TOROIDALTAG);
   }

   template <typename TSimType, template <typename> class TSimTraits> template <typename TVisTraits> void CSCSParaviewGenerator<TSimType, TSimTraits>::writePoloidal()
   {
      // Read data from file
      this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::PoloidalOnly> >();

      // Configure transforms
      this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

      // Transform the fields
      this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

      // Write data to file
      this->mpOutFile->template writeVisualisation<TVisTraits>(CSCSFileDefs::POLOIDALTAG);
   }

   template <typename TSimType, template <typename> class TSimTraits> void CSCSParaviewGenerator<TSimType, TSimTraits>::finalise()
   {
      // Finalise input file IO
      this->mpInFile->finalise();

      // Finish writing to output file
      this->mpOutFile->writeEnd();

      // Finalise output file IO
      this->mpOutFile->finalise();

      // Finalise generator base
      GeneratorBase<TSimType, TSimTraits>::finalise();
   }

}

#endif // CSCSPARAVIEWGENERATOR_HPP
