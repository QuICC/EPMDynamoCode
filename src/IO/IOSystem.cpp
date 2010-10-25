/** \file IOSystem.cpp
 *  \brief Source of the implementation of the IO control system
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/IOSystem.hpp"

// Project includes
//

namespace EPMDynamo {

   IOSystem::IOSystem()
      : mpCfgFile(new ParametersFile("EkQRaRo"))
   {
      // Initialise
      this->init();
   }

   void IOSystem::init()
   {
      // Initialise config file
      this->mpCfgFile->init();
      // Read configuration data from config file
      this->mpCfgFile->read();
      // Close the config file
      this->mpCfgFile->finalise();

      // Initialise the StdMessage ouput
      SmartStdMessage pStd(new StdMessage(mpCfgFile->stdMessageName()));
      this->mpStdMessage = pStd;
      this->mpStdMessage->init();

      // Print the run information
      this->mpCfgFile->printInfo();
   }

   const ArrayI& IOSystem::aTrunc() const
   {
      return this->mpCfgFile->aTrunc();
   }

   const Array& IOSystem::aTStep() const
   {
      return this->mpCfgFile->aTStep();
   }

   const ArrayI& IOSystem::aRunI() const
   {
      return this->mpCfgFile->aRunI();
   }

   const Array& IOSystem::aRun() const
   {
      return this->mpCfgFile->aRun();
   }

   const Array& IOSystem::aEquation() const
   {
      return this->mpCfgFile->aEquation();
   }

   void IOSystem::useInitialState(SmartStateFileReaderBase pInState, TimestepParameters& tsParams)
   {
      // Initialise IO
      pInState->init();

      // Read in data
      pInState->read();

      // Handle time and timestep input
      if(this->aTStep()(0) == -1.0)
      {
         tsParams.setTime(pInState->time());
      }
      if(this->aTStep()(1) == -1.0)
      {
         tsParams.updateTimestep(pInState->timestep());
      }

      // Finalise IO
      pInState->finalise();
   }

   void IOSystem::useImposedField(SmartImposedFieldReaderBase pIField)
   {
      // Initialise IO
      pIField->init();

      // Read file
      pIField->read();

      // Finalise IO
      pIField->finalise();
   }

   void IOSystem::useSource(SmartSourceReaderBase pSource)
   {
      // Initialise IO
      pSource->init();

      // Read file
      pSource->read();

      // Finalise IO
      pSource->finalise();
   }

   void IOSystem::addASCIIWriter(SmartASCIIWriter pWriter)
   {
      this->mASCIIWriters.push_back(pWriter);
   }

   void IOSystem::addHDF5Writer(SmartHDF5Writer pWriter)
   {
      this->mHDF5Writers.push_back(pWriter);
   }

   void IOSystem::initWriters()
   {
      for(unsigned int i=0; i < this->mASCIIWriters.size(); ++i)
      {
         this->mASCIIWriters.at(i)->init();
      }

      for(unsigned int i=0; i < this->mHDF5Writers.size(); ++i)
      {
         this->mHDF5Writers.at(i)->init();
      }
   }

   void IOSystem::finaliseWriters()
   {
      for(unsigned int i=0; i < this->mASCIIWriters.size(); ++i)
      {
         this->mASCIIWriters.at(i)->finalise();
      }

      for(unsigned int i=0; i < this->mHDF5Writers.size(); ++i)
      {
         this->mHDF5Writers.at(i)->finalise();
      }

      this->mpStdMessage->finalise();
   }

   void IOSystem::writeASCII()
   {
      for(unsigned int i=0; i < this->mASCIIWriters.size(); ++i)
      {
         this->mASCIIWriters.at(i)->write();
      }
   }

   void IOSystem::writeHDF5()
   {
      for(unsigned int i=0; i < this->mHDF5Writers.size(); ++i)
      {
         this->mHDF5Writers.at(i)->write();
      }
   }
}
