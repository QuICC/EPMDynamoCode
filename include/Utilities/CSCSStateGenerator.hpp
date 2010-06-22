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

   template <typename TSimType, template <typename> class TGenTraits> void CSCSStateGenerator<TSimType, TGenTraits>::setState(typename CSCSFieldTraits<TSimType, TGenTraits>::CodType &codC)
   {
      // Set state file to read Spectral values from
      EPMSHARED_PTR<CodStateFileReader<TSimType, CSCSFieldTraits> > pInState(new CodStateFileReader<TSimType, CSCSFieldTraits>(codC,  "4242"));

      // Read Spectral values
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create state file for codensity
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<CodCSCSFileWriter<TSimType, CSCSFieldTraits> >  pOutState(new CodCSCSFileWriter<TSimType, CSCSFieldTraits>("FSOuter", codC, this->mSimControl.tsParams()));

      // Set state file as output file and initialise system
      initStateFile(pOutState);
   }

   template <typename TSimType, template <typename> class TGenTraits> void CSCSStateGenerator<TSimType, TGenTraits>::setState(typename CSCSFieldTraits<TSimType, TGenTraits>::MagType  &magB)
   {
      // Set state file to read Spectral values from
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<MagStateFileReader<TSimType, CSCSFieldTraits> > pInState(new MagStateFileReader<TSimType, CSCSFieldTraits>(magB,  "4242"));

      // Read Spectral values
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create state file for codensity
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<MagCSCSFileWriter<TSimType, CSCSFieldTraits> >  pOutState(new MagCSCSFileWriter<TSimType, CSCSFieldTraits>("FSOuter", magB, this->mSimControl.tsParams()));

      // Set state file as output file and initialise system
      initStateFile(pOutState);
   }

   template <typename TSimType, template <typename> class TGenTraits> void CSCSStateGenerator<TSimType, TGenTraits>::setState(typename CSCSFieldTraits<TSimType, TGenTraits>::VelType  &velV)
   {
      // Set state file to read Spectral values from
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<VelStateFileReader<TSimType, CSCSFieldTraits> > pInState(new VelStateFileReader<TSimType, CSCSFieldTraits>(velV,  "4242"));

      // Read Spectral values
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create state file for codensity
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<VelCSCSFileWriter<TSimType, CSCSFieldTraits> >  pOutState(new VelCSCSFileWriter<TSimType, CSCSFieldTraits>("FSOuter", velV, this->mSimControl.tsParams()));

      // Set state file as output file and initialise system
      initStateFile(pOutState);
   }

   template <typename TSimType, template <typename> class TGenTraits> void CSCSStateGenerator<TSimType, TGenTraits>::setState(typename CSCSFieldTraits<TSimType, TGenTraits>::CodType &codC, typename CSCSFieldTraits<TSimType, TGenTraits>::VelType &velV)
   {
      // Set state file to read Spectral values from
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<ThermConvStateFileReader<TSimType, CSCSFieldTraits> > pInState(new ThermConvStateFileReader<TSimType, CSCSFieldTraits>(codC, velV,  "4242"));

      // Read Spectral values
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create state file for codensity
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<ThermConvCSCSFileWriter<TSimType, CSCSFieldTraits> >  pOutState(new ThermConvCSCSFileWriter<TSimType, CSCSFieldTraits>("FSOuter", codC, velV, this->mSimControl.tsParams()));

      // Set state file as output file and initialise system
      initStateFile(pOutState);
   }

   template <typename TSimType, template <typename> class TGenTraits> void CSCSStateGenerator<TSimType, TGenTraits>::setState(typename CSCSFieldTraits<TSimType, TGenTraits>::MagType  &magB, typename CSCSFieldTraits<TSimType, TGenTraits>::VelType &velV)
   {
      // Set state file to read Spectral values from
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<MagConvStateFileReader<TSimType, CSCSFieldTraits> > pInState(new MagConvStateFileReader<TSimType, CSCSFieldTraits>(magB, velV,  "4242"));

      // Read Spectral values
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create state file for codensity
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<MagConvCSCSFileWriter<TSimType, CSCSFieldTraits> >  pOutState(new MagConvCSCSFileWriter<TSimType, CSCSFieldTraits>("FSOuter", magB, velV, this->mSimControl.tsParams()));

      // Set state file as output file and initialise system
      initStateFile(pOutState);
   }

   template <typename TSimType, template <typename> class TGenTraits> void CSCSStateGenerator<TSimType, TGenTraits>::setState(typename CSCSFieldTraits<TSimType, TGenTraits>::CodType &codC, typename CSCSFieldTraits<TSimType, TGenTraits>::MagType &magB, typename CSCSFieldTraits<TSimType, TGenTraits>::VelType &velV)
   {
      // Set state file to read Spectral values from
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<StateFileReader<TSimType, CSCSFieldTraits> > pInState(new StateFileReader<TSimType, CSCSFieldTraits>(codC, magB, velV,  "4242"));

      // Read Spectral values
      this->mIOSys.useInitialState(pInState, this->mSimControl.tsParams());

      // Create state file for codensity
      EPMDYNAMO_SHAREDPTRNS::shared_ptr<CSCSFileWriter<TSimType, CSCSFieldTraits> >  pOutState(new CSCSFileWriter<TSimType, CSCSFieldTraits>("FSOuter", codC, magB, velV, this->mSimControl.tsParams()));

      // Set state file as output file and initialise system
      initStateFile(pOutState);
   }

}

#endif // CSCSSTATEGENERATOR_HPP
