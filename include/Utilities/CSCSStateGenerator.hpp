/** \file CSCSStateGenerator.hpp
 *  \brief Implementation of a simple CSCS visualisation file format generator
 */

#ifndef CSCSSTATEGENERATOR_HPP
#define CSCSSTATEGENERATOR_HPP

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
#include "IO/HDF5/Visualisation/CSCSFileWriter.hpp"
#include "IO/HDF5/State/StateFileReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a simple CSCS visualisation file format generator
    *
    * \tparam TSimType Type of the simulation
    * \tparam TGenTraits Traits for the generator
    */
   template <typename TSimType, template <typename> class TGenTraits> class CSCSStateGenerator: public GeneratorBase<TSimType, TGenTraits>
   {
      public:
         /**
          * @brief Constructor
          */
         CSCSStateGenerator();

         /**
          * @brief Simple empty destructor
          */
         virtual ~CSCSStateGenerator() {};

         /**
          * @brief Initialise the output file
          */
         virtual void initOutput(std::string name = "CSCSParaView");

      protected:

      private:
   };

   template <typename TSimType, template <typename> class TGenTraits> CSCSStateGenerator<TSimType, TGenTraits>::CSCSStateGenerator()
      : GeneratorBase<TSimType, TGenTraits>()
   {
   }

   template <typename TSimType, template <typename> class TGenTraits> void CSCSStateGenerator<TSimType, TGenTraits>::initOutput(std::string name)
   {
      EPMSHARED_PTR<StateFileReader<TSimType, TGenTraits> > pInFile;
      EPMSHARED_PTR<StateFileWriter<TSimType, TGenTraits> >  pOutFile;

      // Create file with all fields
      if(TGenTraits<TSimType>::NeedCodensity && TGenTraits<TSimType>::NeedMagnetic && TGenTraits<TSimType>::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimType, TGenTraits>(this->codC(), this->magB(), this->velV(), "4242"));

         pOutFile.reset(new CSCSStateFileWriter<TSimType, TGenTraits>("FSOuter", this->codC(), this->magB(), this->velV(), this->mSimControl.tsParams()));

      // Create file with Codensity and Velocity fields
      } else if(TGenTraits<TSimType>::NeedCodensity && TGenTraits<TSimType>::NeedMagnetic)
      {
         pInFile.reset(new StateFileReader<TSimType, TGenTraits>(this->codC(), this->velV(), "4242"));

         pOutFile.reset(new CSCSStateFileWriter<TSimType, TGenTraits>("FSOuter", this->codC(), this->velV(), this->mSimControl.tsParams()));

      // Create file with Magnetic and Velocity fields
      } else if(TGenTraits<TSimType>::NeedMagnetic && TGenTraits<TSimType>::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimType, TGenTraits>(this->magB(), this->velV(), "4242"));

         pOutFile.reset(new CSCSStateFileWriter<TSimType, TGenTraits>("FSOuter", this->magB(), this->velV(), this->mSimControl.tsParams()));

      // Create file with only codensity field
      } else if(TGenTraits<TSimType>::NeedCodensity)
      {
         pInFile.reset(new StateFileReader<TSimType, TGenTraits>(this->codC(), "4242"));

         pOutFile.reset(new CSCSStateFileWriter<TSimType, TGenTraits>("FSOuter", this->codC(), this->mSimControl.tsParams()));

      // Create file with only magnetic field
      } else if(TGenTraits<TSimType>::NeedMagnetic)
      {
         pInFile.reset(new StateFileReader<TSimType, TGenTraits>(this->magB(), "4242"));

         pOutFile.reset(new CSCSStateFileWriter<TSimType, TGenTraits>("FSOuter", this->magB(), this->mSimControl.tsParams()));

      // Create file with only velocity field
      } else if(TGenTraits<TSimType>::NeedVelocity)
      {
         pInFile.reset(new StateFileReader<TSimType, TGenTraits>(this->velV(), "4242"));

         pOutFile.reset(new CSCSStateFileWriter<TSimType, TGenTraits>("FSOuter", this->velV(), this->mSimControl.tsParams()));
      }

      // Read Spectral values
      this->mIOSys.useInitialState(pInFile, this->mSimControl.tsParams());

      // Change the base name
      pOutFile->changeBasename(name);

      // Set state file as output file and initialise system
      this->initOutputFile(pOutFile);
   }

}

#endif // CSCSSTATEGENERATOR_HPP
