/** \file StdMessage.cpp
 *  \brief Source of the implementation of the std message output
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "IO/ASCII/StdMessage.hpp"

// Project includes
//
#include "IO/ASCII/StdMessageDefs.hpp"

namespace EPMDynamo {

   StdMessage::StdMessage(std::string name)
      : ASCIIEWriter(name + StdMessageDefs::BASENAME, StdMessageDefs::EXTENSION, StdMessageDefs::HEADER, StdMessageDefs::VERSION)
   {
   }

   void StdMessage::init()
   {
      ASCIIEWriter::init();

      // Backup std::cout buffer
      this->mpCoutBuffer = std::cout.rdbuf();

      // Redirect std::cout to this file
      std::cout.rdbuf(this->mFile.rdbuf());
   }

   void StdMessage::write()
   {
   }

   void StdMessage::finalise()
   {
      // Put std::cout buffer back in place
      std::cout.rdbuf(this->mpCoutBuffer);

      // close file etc
      ASCIIEWriter::finalise();
   }
}
