/** \file ImposedFieldGenerator.hpp
 *  \brief Implementation of a imposed field file generator
 */

#ifndef IMPOSEDFIELDGENERATOR_HPP
#define IMPOSEDFIELDGENERATOR_HPP

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
#include "IO/HDF5/Imposed/ImposedFieldWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a few methods for easier implementation of an initial state generator
    *
    * \tparam TSimType Type of the simulation
    * \tparam TSimTraits Traits for the generator
    */
   template <typename TSimType, template <typename> class TSimTraits> class ImposedFieldGenerator: public GeneratorBase<TSimType, TSimTraits>
   {
      public:
         /// Typedef for a smart state file
         typedef  EPMSHARED_PTR<ImposedFieldWriter<TSimType, TSimTraits> >   SmartImposedWriter; 

         /**
          * @brief Constructor
          */
         ImposedFieldGenerator();

         /**
          * @brief Simple empty destructor
          */
         virtual ~ImposedFieldGenerator() {};

         /**
          * @brief Initialise the output file
          *
          * @param name Basename of the output file
          */
         virtual void setupOutput(std::string name = "Imposed");

         /**
          * @brief Write the imposed file
          */
         void writeImposedFile();

         /**
          * @brief Finalise the output imposed file
          */
         void finalise();

      protected:
         /**
          * @brief Smart pointer for the imposed field file
          */
         SmartImposedWriter   mpOutFile;

         /**
          * @brief Initialise the output imposed field file
          *
          * @param pFile Output file smart pointer
          */
         void initImposedFile(SmartImposedWriter pFile);

      private:
   };

   template <typename TSimType, template <typename> class TSimTraits> ImposedFieldGenerator<TSimType, TSimTraits>::ImposedFieldGenerator()
      : GeneratorBase<TSimType, TSimTraits>()
   {
   }

   template <typename TSimType, template <typename> class TSimTraits> void ImposedFieldGenerator<TSimType, TSimTraits>::setupOutput(std::string name)
   {
      SmartImposedWriter  pOutFile;

      // Create state with all fields
      if(TSimTraits<TSimType>::NeedCodensity && TSimTraits<TSimType>::NeedMagnetic && TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new ImposedFieldWriter<TSimType, TSimTraits>(this->codC(), this->magB(), this->velV()));

      // Create state with Codensity and Velocity fields
      } else if(TSimTraits<TSimType>::NeedCodensity && TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new ImposedFieldWriter<TSimType, TSimTraits>(this->codC(), this->velV()));

      // Create state with Magnetic and Velocity fields
      } else if(TSimTraits<TSimType>::NeedMagnetic && TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new ImposedFieldWriter<TSimType, TSimTraits>(this->magB(), this->velV()));

      // Create state with only codensity field
      } else if(TSimTraits<TSimType>::NeedCodensity)
      {
         pOutFile.reset(new ImposedFieldWriter<TSimType, TSimTraits>(this->codC()));

      // Create state with only magnetic field
      } else if(TSimTraits<TSimType>::NeedMagnetic)
      {
         pOutFile.reset(new ImposedFieldWriter<TSimType, TSimTraits>(this->magB()));

      // Create state with only velocity field
      } else if(TSimTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new ImposedFieldWriter<TSimType, TSimTraits>(this->velV()));
      }

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initImposedFile(pOutFile);
   }

   template <typename TSimType, template <typename> class TSimTraits> void ImposedFieldGenerator<TSimType, TSimTraits>::initImposedFile(typename ImposedFieldGenerator<TSimType, TSimTraits>::SmartImposedWriter pFile)
   {
      // Set output file
      this->mpOutFile = pFile;
   }

   template <typename TSimType, template <typename> class TSimTraits> void ImposedFieldGenerator<TSimType, TSimTraits>::writeImposedFile()
   {
      this->mpOutFile->write();
   }

   template <typename TSimType, template <typename> class TSimTraits> void ImposedFieldGenerator<TSimType, TSimTraits>::finalise()
   {
      // Finalise output file
      this->mpOutFile->finalise();

      // Finalise generator
      GeneratorBase<TSimType, TSimTraits>::finalise();
   }

}

#endif // IMPOSEDFIELDGENERATOR_HPP
