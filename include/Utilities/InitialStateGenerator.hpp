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
    * \tparam TGenTraits Traits for the generator
    */
   template <typename TSimType, template <typename> class TGenTraits> class InitialStateGenerator: public GeneratorBase<TSimType, TGenTraits>
   {
      public:
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
          */
         virtual void initOutput(std::string name = "Initial");

      protected:

      private:
   };

   template <typename TSimType, template <typename> class TGenTraits> InitialStateGenerator<TSimType, TGenTraits>::InitialStateGenerator()
      : GeneratorBase<TSimType, TGenTraits>()
   {
   }

   template <typename TSimType, template <typename> class TGenTraits> void InitialStateGenerator<TSimType, TGenTraits>::initOutput(std::string name)
   {
      EPMSHARED_PTR<StateFileWriter<TSimType, TGenTraits> >  pOutFile;

      // Create state with all fields
      if(TGenTraits<TSimType>::NeedCodensity && TGenTraits<TSimType>::NeedMagnetic && TGenTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TGenTraits>(this->codC(), this->magB(), this->velV(), this->mEqParams, this->mSimControl.tsParams()));

      // Create state with Codensity and Velocity fields
      } else if(TGenTraits<TSimType>::NeedCodensity && TGenTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TGenTraits>(this->codC(), this->velV(), this->mEqParams, this->mSimControl.tsParams()));

      // Create state with Magnetic and Velocity fields
      } else if(TGenTraits<TSimType>::NeedMagnetic && TGenTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TGenTraits>(this->magB(), this->velV(), this->mEqParams, this->mSimControl.tsParams()));

      // Create state with only codensity field
      } else if(TGenTraits<TSimType>::NeedCodensity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TGenTraits>(this->codC(), this->mEqParams, this->mSimControl.tsParams()));

      // Create state with only magnetic field
      } else if(TGenTraits<TSimType>::NeedMagnetic)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TGenTraits>(this->magB(), this->mEqParams, this->mSimControl.tsParams()));

      // Create state with only velocity field
      } else if(TGenTraits<TSimType>::NeedVelocity)
      {
         pOutFile.reset(new StateFileWriter<TSimType, TGenTraits>(this->velV(), this->mEqParams, this->mSimControl.tsParams()));
      }

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initOutputFile(pOutFile);
   }

}

#endif // INITIALSTATEGENERATOR_HPP
