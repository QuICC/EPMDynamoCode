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


#include <fstream>

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
          *
          * @param type Type of the additional
          */
         template <typename TAddVisTraits> bool prepareAdditional(const int type);

         /**
          * @brief Write total field values
          *
          * @param type Type of the additional
          */
         template <typename TAddVisTraits> void finaliseAdditional(const int type);

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

   template <typename TSimTraits> template <typename TAddVisTraits> bool CSCSParaviewGenerator<TSimTraits>::prepareAdditional(const int type)
   {
      bool needSave = false;

      if(type == 1 && this->mpVelV != NULL && TAddVisTraits::VisInertial)
      {
         Array amplitude;

         amplitude = this->velV().rOc().rPerturbation().computeXYSolidProjection(this->mTransform.radBasis());

         // Output solid body projection values before operation
         std::cerr << "Solid body projection before: " << amplitude.transpose() << std::endl;

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

         // Output solid body projection values after removal
         std::cerr << "Solid body projection after: "  << amplitude.transpose() << std::endl;

         needSave = true;
      }

      if(type == 2 && this->mpVelV != NULL && TAddVisTraits::VisVelM0)
      {
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
                  if(m_ != 0)
                  {
                    this->velV().rOc().rPerturbation().rTor().rLShell(l).col(m).setConstant(0.0);
                    this->velV().rOc().rPerturbation().rPol().rLShell(l).col(m).setConstant(0.0);
                  }
               }
            }
         }

         needSave = true;
      }

      if(type == 2 && this->mpMagB != NULL && TAddVisTraits::VisMagM0)
      {
         // Get truncation information
         SmartTruncation pTrunc = this->magB().oc().trunc();

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
                  if(m_ != 0)
                  {
                    this->magB().rOc().rPerturbation().rTor().rLShell(l).col(m).setConstant(0.0);
                    this->magB().rOc().rPerturbation().rPol().rLShell(l).col(m).setConstant(0.0);
                  }
               }
            }
         }

         needSave = true;
      }

      return needSave;
   }

   template <typename TSimTraits> template <typename TAddVisTraits> void CSCSParaviewGenerator<TSimTraits>::finaliseAdditional(const int type)
   {
      if(type == 1 && this->mpVelV != NULL && TAddVisTraits::VisInertial)
      {
      }

      if(type == 2 && this->mpVelV != NULL && TAddVisTraits::VisVelM0)
      {
      }

      if(type == 2 && this->mpVelV != NULL && TAddVisTraits::VisASCIIVelM0)
      {
         SmartTruncation pTrunc = this->velV().oc().trunc();
         std::ofstream file;

         file.open("VelM0_gridR.dat");
         file << pTrunc->sim()->rad()->radGrid();
         file.close();
         file.open("VelM0_gridTh.dat");
         file << pTrunc->sim()->hoz()->thGrid();
         file.close();

         file.open("VelM0_radial.dat");

         for(int r = 0; r < pTrunc->local()->rtp()->nR(); r++)
         {
            file << this->velV().oc().rtp().r().shell(r).row(0) << std::endl;
         }

         file.close();

         file.open("VelM0_theta.dat");

         for(int r = 0; r < pTrunc->local()->rtp()->nR(); r++)
         {
            file << this->velV().oc().rtp().theta().shell(r).row(0) << std::endl;
         }

         file.close();

         file.open("VelM0_phi.dat");

         for(int r = 0; r < pTrunc->local()->rtp()->nR(); r++)
         {
            file << this->velV().oc().rtp().phi().shell(r).row(0) << std::endl;
         }

         file.close();
      }

      if(type == 2 && this->mpVelV != NULL && TAddVisTraits::VisASCIIVel)
      {
         SmartTruncation pTrunc = this->velV().oc().trunc();
         std::ofstream file;

         file.open("Vel_gridR.dat");
         file << pTrunc->sim()->rad()->radGrid();
         file.close();
         file.open("Vel_gridTh.dat");
         file << pTrunc->sim()->hoz()->thGrid();
         file.close();
         file.open("Vel_gridPh.dat");
         file << pTrunc->sim()->hoz()->phGrid();
         file.close();

         file.open("Vel_radial.dat");

         for(int r = 0; r < pTrunc->local()->rtp()->nR(); r++)
         {
            file << this->velV().oc().rtp().r().shell(r) << std::endl;
         }

         file.close();

         file.open("Vel_theta.dat");

         for(int r = 0; r < pTrunc->local()->rtp()->nR(); r++)
         {
            file << this->velV().oc().rtp().theta().shell(r) << std::endl;
         }

         file.close();

         file.open("Vel_phi.dat");

         for(int r = 0; r < pTrunc->local()->rtp()->nR(); r++)
         {
            file << this->velV().oc().rtp().phi().shell(r) << std::endl;
         }

         file.close();
      }

      if(type == 2 && this->mpMagB != NULL && TAddVisTraits::VisMagM0)
      {
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


      bool hasAdditional;
      for(int i=0; i < 3; i++)
      {
         // initialise flag
         hasAdditional = false;

         // Read data from file
         this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::FullField> >();

         hasAdditional = this->template prepareAdditional<TAddVisTraits>(i);

         // Configure transforms
         this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

         // Transform the fields
         this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

         this->template finaliseAdditional<TAddVisTraits>(i);

         // Write data to file
         this->mpOutFile->template writeAdditional<TAddVisTraits>(i);
      }
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


      bool hasAdditional;
      for(int i=0; i < 3; i++)
      {
         // initialise flag
         hasAdditional = false;

         // Read data from file
         this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::ToroidalOnly> >();

         hasAdditional = this->template prepareAdditional<TAddVisTraits>(i);

         // Configure transforms
         this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

         // Transform the fields
         this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

         this->template finaliseAdditional<TAddVisTraits>(i);

         // Write data to file
         this->mpOutFile->template writeAdditional<TAddVisTraits>(i,CSCSFileDefs::TOROIDALTAG);
      }
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


      bool hasAdditional;
      for(int i=0; i < 3; i++)
      {
         // initialise flag
         hasAdditional = false;

         // Read data from file
         this->mpInFile->template readPartial<VisStateFilterTraits<TVisTraits, StateFileDefs::PoloidalOnly> >();

         hasAdditional = this->template prepareAdditional<TAddVisTraits>(i);

         // Configure transforms
         this->template configureTransform<VisGeneratorTraits<TVisTraits> >();

         // Transform the fields
         this->template transformSpectral<VisGeneratorTraits<TVisTraits> >();

         this->template finaliseAdditional<TAddVisTraits>(i);

         // Write data to file
         this->mpOutFile->template writeAdditional<TAddVisTraits>(i,CSCSFileDefs::POLOIDALTAG);
      }
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
