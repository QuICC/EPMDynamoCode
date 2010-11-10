/** \file IOBaseSystem.cpp
 *  \brief Source of the base properties of the IO control system
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/IOBaseSystem.hpp"

// Project includes
//

namespace EPMDynamo {

   IOBaseSystem::IOBaseSystem()
      : mpCfgFile(new ParametersFile("EkQRaRo"))
   {
      // Initialise
      this->initBase();
   }

   void IOBaseSystem::initBase()
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

   void IOBaseSystem::finaliseBase()
   {
      this->mpStdMessage->finalise();
   }
}
