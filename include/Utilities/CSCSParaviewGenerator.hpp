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
    * \tparam TSimTraits Traits for the generator
    */
   template <typename TSimTraits> class CSCSParaviewGenerator: public GeneratorBase<TSimTraits>
   {
      public:
         /// Typedef for the state file reader
         typedef EPMSHARED_PTR<StateFileReader<TSimTraits> >  SmartStateReader; 

         /// Typedef for the CSCS file writer
         typedef EPMSHARED_PTR<CSCSFileWriter<TSimTraits> >   SmartCSCSWriter;

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
         template <typename TVisTraits, typename TAddVisTraits> void writeTotal();

         /**
          * @brief Write toroidal component values
          */
         template <typename TVisTraits, typename TAddVisTraits> void writeToroidal();

         /**
          * @brief Write poloidal component values
          */
         template <typename TVisTraits, typename TAddVisTraits> void writePoloidal();

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

         /**
          * @brief Write total field values
          */
         template <typename TAddVisTraits> void prepareAdditional();

      private:
   };

   template <typename TSimTraits> CSCSParaviewGenerator<TSimTraits>::CSCSParaviewGenerator()
      : GeneratorBase<TSimTraits>()
   {
   }

   template <typename TSimTraits> void CSCSParaviewGenerator<TSimTraits>::setupOutput(std::string name)
   {
      SmartStateReader  pInFile;
      SmartCSCSWriter  pOutFile;

      // Create file with all fields
      if(TSimTraits::NeedCodensity && TSimTraits::NeedMagnetic && TSimTraits::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->codC(), this->magB(), this->velV(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimTraits>("FSOuter", this->codC(), this->magB(), this->velV(), this->mTSParams));

      // Create file with Codensity and Velocity fields
      } else if(TSimTraits::NeedCodensity && TSimTraits::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->codC(), this->velV(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimTraits>("FSOuter", this->codC(), this->velV(), this->mTSParams));

      // Create file with Magnetic and Velocity fields
      } else if(TSimTraits::NeedMagnetic && TSimTraits::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->magB(), this->velV(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimTraits>("FSOuter", this->magB(), this->velV(), this->mTSParams));

      // Create file with only codensity field
      } else if(TSimTraits::NeedCodensity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->codC(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimTraits>("FSOuter", this->codC(), this->mTSParams));

      // Create file with only magnetic field
      } else if(TSimTraits::NeedMagnetic)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->magB(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimTraits>("FSOuter", this->magB(), this->mTSParams));

      // Create file with only velocity field
      } else if(TSimTraits::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimTraits>(this->velV(), "4CSCS"));

         pOutFile.reset(new CSCSFileWriter<TSimTraits>("FSOuter", this->velV(), this->mTSParams));
      }

      // Init the input file
      this->initInputFile(pInFile);

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initCSCSFile(pOutFile);
   }

   template <typename TSimTraits> void CSCSParaviewGenerator<TSimTraits>::initInputFile(typename CSCSParaviewGenerator<TSimTraits>::SmartStateReader   pInFile)
   {
      // Set the state file localy
      this->mpInFile = pInFile;

      // Initialise the input file
      this->mpInFile->init();

      // Read setup information from file
      this->mpInFile->readSetup();
   }

   template <typename TSimTraits> void CSCSParaviewGenerator<TSimTraits>::initCSCSFile(typename CSCSParaviewGenerator<TSimTraits>::SmartCSCSWriter   pFile)
   {
      // Set output file
      this->mpOutFile = pFile;

      // Set output file
      this->mpOutFile->init();

      // Write start information to file
      this->mpOutFile->writeStart();
   }

   template <typename TSimTraits> template <typename TAddVisTraits> void CSCSParaviewGenerator<TSimTraits>::prepareAdditional()
   {
      if(this->mpVelV != NULL && TAddVisTraits::VisInertial)
      {
         Array amplitude;

         amplitude = this->velV().rOc().rPerturbation().computeXYSolidProjection(this->mTransform.radBasis());

         std::cerr << amplitude.transpose() << std::endl;

         // Get truncation information
         SmartTruncation pTrunc = this->velV().oc().trunc();

         // Get the stored ls
         ArrayI ls = pTrunc->local()->spec()->lArray();
         ArrayI ms;
         int l_ = -1;
         int m_ = -1;
         // Loop over all stored l's
         for(int l = 0; l < ls.size(); l++)
         {
            l_ = ls(l);
            if(l_ == 1)
            {
               ms = pTrunc->local()->spec()->mArray(l);
               for(int m =0; m < ms.size(); m++)
               {
                  m_ = ms(m);
                  if(m_ == 1)
                  {
                    this->velV().rOc().rPerturbation().rTor().rLShell(l).col(m)(0).real() -= amplitude(0);
                    this->velV().rOc().rPerturbation().rTor().rLShell(l).col(m)(0).imag() -= amplitude(1);
                  }
               }
            }
         }

         amplitude = this->velV().rOc().rPerturbation().computeXYSolidProjection(this->mTransform.radBasis());

         std::cerr << amplitude.transpose() << std::endl;
      }
   }

   template <typename TSimTraits> template <typename TVisTraits, typename TAddVisTraits> void CSCSParaviewGenerator<TSimTraits>::writeTotal()
   {
      // Read data from file
      this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::FullField> >();

      // Configure transforms
      this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

      // Transform the fields
      this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

      // Write data to file
      this->mpOutFile->template writeVisualisation<TVisTraits>();



      // Read data from file
      this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::FullField> >();

      this->template prepareAdditional<TAddVisTraits>();

      // Configure transforms
      this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

      // Transform the fields
      this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

      // Write data to file
      this->mpOutFile->template writeAdditional<TAddVisTraits>();
   }

   template <typename TSimTraits> template <typename TVisTraits, typename TAddVisTraits> void CSCSParaviewGenerator<TSimTraits>::writeToroidal()
   {
      // Read data from file
      this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::ToroidalOnly> >();

      // Configure transforms
      this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

      // Transform the fields
      this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

      // Write data to file
      this->mpOutFile->template writeVisualisation<TVisTraits>(CSCSFileDefs::TOROIDALTAG);


      // Read data from file
      this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::ToroidalOnly> >();

      this->template prepareAdditional<TAddVisTraits>();

      // Configure transforms
      this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

      // Transform the fields
      this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

      // Write data to file
      this->mpOutFile->template writeAdditional<TAddVisTraits>(CSCSFileDefs::TOROIDALTAG);
   }

   template <typename TSimTraits> template <typename TVisTraits, typename TAddVisTraits> void CSCSParaviewGenerator<TSimTraits>::writePoloidal()
   {
      // Read data from file
      this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::PoloidalOnly> >();

      // Configure transforms
      this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

      // Transform the fields
      this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

      // Write data to file
      this->mpOutFile->template writeVisualisation<TVisTraits>(CSCSFileDefs::POLOIDALTAG);



      // Read data from file
      this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::PoloidalOnly> >();

      this->template prepareAdditional<TAddVisTraits>();

      // Configure transforms
      this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

      // Transform the fields
      this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

      // Write data to file
      this->mpOutFile->template writeAdditional<TAddVisTraits>(CSCSFileDefs::POLOIDALTAG);
   }

   template <typename TSimTraits> void CSCSParaviewGenerator<TSimTraits>::finalise()
   {
      // Finalise input file IO
      this->mpInFile->finalise();

      // Finish writing to output file
      this->mpOutFile->writeEnd();

      // Finalise output file IO
      this->mpOutFile->finalise();

      // Finalise generator base
      GeneratorBase<TSimTraits>::finalise();
   }

}

#endif // CSCSPARAVIEWGENERATOR_HPP
