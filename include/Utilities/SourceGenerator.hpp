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
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits for the generator
    */
   template <typename TSimType, template <typename> class TSimTraits> class SourceGenerator: public GeneratorBase<TSimType, TSimTraits>
   {
      public:
         /// Typedef for a smart state file
         typedef  EPMSHARED_PTR<SourceFileWriter<TSimType, TSimTraits> >   SmartSourceWriter; 

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

   template <typename TSimType, template <typename> class TSimTraits> SourceGenerator<TSimType, TSimTraits>::SourceGenerator()
      : GeneratorBase<TSimType, TSimTraits>()
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void SourceGenerator<TSimType, TSimTraits>::setupOutput(std::string name)
   {
      SmartSourceWriter pOutFile;

      // Create state with all fields
      if(TSimTraits<TSimType>::NeedCodensity && TSimTraits<TSimType>::NeedMagnetic && TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new SourceFileWriter<TSimType, TSimTraits>(this->codC(), this->magB(), this->velV()));

      // Create state with Codensity and Velocity fields
      } else if(TSimTraits<TSimType>::NeedCodensity && TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new SourceFileWriter<TSimType, TSimTraits>(this->codC(), this->velV()));

      // Create state with Magnetic and Velocity fields
      } else if(TSimTraits<TSimType>::NeedMagnetic && TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new SourceFileWriter<TSimType, TSimTraits>(this->magB(), this->velV()));

      // Create state with only codensity field
      } else if(TSimTraits<TSimType>::NeedCodensity)
      {
         pOutFile.reset(new SourceFileWriter<TSimType, TSimTraits>(this->codC()));

      // Create state with only magnetic field
      } else if(TSimTraits<TSimType>::NeedMagnetic)
      {
         pOutFile.reset(new SourceFileWriter<TSimType, TSimTraits>(this->magB()));

      // Create state with only velocity field
      } else if(TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new SourceFileWriter<TSimType, TSimTraits>(this->velV()));
      }

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initSourceFile(pOutFile);
   }

   template <typename TSimType, template <typename> class TSimTraits> void SourceGenerator<TSimType, TSimTraits>::initSourceFile(typename SourceGenerator<TSimType, TSimTraits>::SmartSourceWriter pFile)
   {
      // Set output file
      this->mpOutFile = pFile;
   }

   template <typename TSimType, template <typename> class TSimTraits> void SourceGenerator<TSimType, TSimTraits>::writeSourceFile()
   {
      this->mpOutFile->write();
   }

   template <typename TSimType, template <typename> class TSimTraits> void SourceGenerator<TSimType, TSimTraits>::finalise()
   {
      // Finalise output file
      this->mpOutFile->finalise();

      // Finalise generator
      GeneratorBase<TSimType, TSimTraits>::finalise();
   }

}

#endif // SOURCEGENERATOR_HPP
