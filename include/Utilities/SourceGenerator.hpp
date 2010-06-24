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
#include "IO/HDF5/Source/SourceWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a few methods for easier implementation of an initial state generator
    *
    * \tparam TSimType Type of the simulation
    * \tparam TGenTraits Traits for the generator
    */
   template <typename TSimType, template <typename> class TGenTraits> class SourceGenerator: public GeneratorBase<TSimType, TGenTraits>
   {
      public:
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
          */
         virtual void initOutput(std::string name = "Source");

      protected:

      private:
   };

   template <typename TSimType, template <typename> class TGenTraits> SourceGenerator<TSimType, TGenTraits>::SourceGenerator()
      : GeneratorBase<TSimType, TGenTraits>()
   {
   }

   template <typename TSimType, template <typename> class TGenTraits> void SourceGenerator<TSimType, TGenTraits>::initOutput(std::string name)
   {
      EPMSHARED_PTR<SourceWriter<TSimType, TGenTraits> >  pOutFile;

      // Create state with all fields
      if(TGenTraits<TSimType>::NeedCodensity && TGenTraits<TSimType>::NeedMagnetic && TGenTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new SourceWriter<TSimType, TGenTraits>(this->codC(), this->magB(), this->velV()));

      // Create state with Codensity and Velocity fields
      } else if(TGenTraits<TSimType>::NeedCodensity && TGenTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new SourceWriter<TSimType, TGenTraits>(this->codC(), this->velV()));

      // Create state with Magnetic and Velocity fields
      } else if(TGenTraits<TSimType>::NeedMagnetic && TGenTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new SourceWriter<TSimType, TGenTraits>(this->magB(), this->velV()));

      // Create state with only codensity field
      } else if(TGenTraits<TSimType>::NeedCodensity)
      {
         pOutFile.reset(new SourceWriter<TSimType, TGenTraits>(this->codC()));

      // Create state with only magnetic field
      } else if(TGenTraits<TSimType>::NeedMagnetic)
      {
         pOutFile.reset(new SourceWriter<TSimType, TGenTraits>(this->magB()));

      // Create state with only velocity field
      } else if(TGenTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new SourceWriter<TSimType, TGenTraits>(this->velV()));
      }

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initOutputFile(pOutFile);
   }

}

#endif // SOURCEGENERATOR_HPP
