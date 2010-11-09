/** \file IOSystemBase.cpp
 *  \brief Source of the base properties of the IO control system
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/IOSystemBase.hpp"

// Project includes
//

namespace EPMDynamo {

   IOSystemBase::IOSystemBase()
      : mpCfgFile(new ParametersFile("EkQRaRo"))
   {
      // Initialise
      this->init();
   }

   void IOSystemBase::init()
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

   void IOSystemBase::finaliseBase()
   {
      this->mpStdMessage->finalise();
   }
}
