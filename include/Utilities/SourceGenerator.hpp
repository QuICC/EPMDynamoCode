/** \file SourceGenerator.hpp
 *  \brief Implementation of a source file generator
 */

#ifndef SOURCEGENERATOR_HPP
#define SOURCEGENERATOR_HPP

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
#include "IO/HDF5/Source/SourceFileWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a few methods for easier implementation of an initial state generator
    *
    * \tparam TSimTraits Traits for the generator
    */
   template <typename TSimTraits> class SourceGenerator: public GeneratorBase<TSimTraits>
   {
      public:
         /// Typedef for a smart state file
         typedef  EPMSHARED_PTR<SourceFileWriter<TSimTraits> >   SmartSourceWriter; 

         /**
          * @brief Constructor
          */
         SourceGenerator();

         /**
          * @brief Simple empty destructor
          */
         virtual ~SourceGenerator() {};

         /**
          * @brief Initialise the output file
          *
          * @param name Basename of the output file
          */
         virtual void setupOutput(std::string name = "Source");

         /**
          * @brief Write the source file
          */
         void writeSourceFile();

         /**
          * @brief Finalise the output source file
          */
         void finalise();

      protected:
         /**
          * @brief Smart pointer for the source file
          */
         SmartSourceWriter   mpOutFile;

         /**
          * @brief Initialise the output source file
          *
          * @param pFile Output file smart pointer
          */
         void initSourceFile(SmartSourceWriter pFile);

      private:
   };

   template <typename TSimTraits> SourceGenerator<TSimTraits>::SourceGenerator()
      : GeneratorBase<TSimTraits>()
   {
   }

   template <typename TSimTraits> void SourceGenerator<TSimTraits>::setupOutput(std::string name)
   {
      SmartSourceWriter pOutFile;

      // Create state with all fields
      if(TSimTraits::NeedCodensity && TSimTraits::NeedMagnetic && TSimTraits::NeedVelocity)
      {
         pOutFile.reset(new SourceFileWriter<TSimTraits>(this->codC(), this->magB(), this->velV()));

      // Create state with Codensity and Velocity fields
      } else if(TSimTraits::NeedCodensity && TSimTraits::NeedVelocity)
      {
         pOutFile.reset(new SourceFileWriter<TSimTraits>(this->codC(), this->velV()));

      // Create state with Magnetic and Velocity fields
      } else if(TSimTraits::NeedMagnetic && TSimTraits::NeedVelocity)
      {
         pOutFile.reset(new SourceFileWriter<TSimTraits>(this->magB(), this->velV()));

      // Create state with only codensity field
      } else if(TSimTraits::NeedCodensity)
      {
         pOutFile.reset(new SourceFileWriter<TSimTraits>(this->codC()));

      // Create state with only magnetic field
      } else if(TSimTraits::NeedMagnetic)
      {
         pOutFile.reset(new SourceFileWriter<TSimTraits>(this->magB()));

      // Create state with only velocity field
      } else if(TSimTraits::NeedVelocity)
      {
         pOutFile.reset(new SourceFileWriter<TSimTraits>(this->velV()));
      }

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initSourceFile(pOutFile);
   }

   template <typename TSimTraits> void SourceGenerator<TSimTraits>::initSourceFile(typename SourceGenerator<TSimTraits>::SmartSourceWriter pFile)
   {
      // Set output file
      this->mpOutFile = pFile;
   }

   template <typename TSimTraits> void SourceGenerator<TSimTraits>::writeSourceFile()
   {
      this->mpOutFile->write();
   }

   template <typename TSimTraits> void SourceGenerator<TSimTraits>::finalise()
   {
      // Finalise output file
      this->mpOutFile->finalise();

      // Finalise generator
      GeneratorBase<TSimTraits>::finalise();
   }

}

#endif // SOURCEGENERATOR_HPP
