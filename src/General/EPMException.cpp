/** \file EPMException.cpp
 *  \brief Definitions of the EPMException methods.
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "General/EPMException.hpp"

// Project includes
//

namespace EPMDynamo {

   EPMException::EPMException(std::string location, std::string message)
      : mLocation(location), mMessage(message)   
   {
   }

   void EPMException::printStdMessage() const
   {
      std::cerr << "***************************************************" << std::endl;
      std::cerr << "EPMException has been caught" << std::endl;
      std::cerr << "    at: " << this->location() << std::endl;
      std::cerr << "    msg: " << this->error() << std::endl;
      std::cerr << "***************************************************" << std::endl;
   }

}
