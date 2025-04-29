/** \file ASCIIRWriter.cpp
 *  \brief Source of the replacing ASCII writer
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/ASCII/ASCIIRWriter.hpp"

// Project includes
//

namespace EPMDynamo {

   ASCIIRWriter::ASCIIRWriter(std::string name, std::string ext, std::string header, std::string version)
      : ASCIIWriter(name, ext, header, version)
   {
   }

   void ASCIIRWriter::init()
   {
   }

   void ASCIIRWriter::finalise()
   {
   }

   void ASCIIRWriter::preWrite()
   {
      if(this->doesIO())
      {
         // Create file
         this->open();

         // Add header information
         this->mFile << FILE_HEADER + this->mHeader << std::endl;
         this->mFile << FILE_VERSION + " " + this->mVersion << std::endl;
      }
   }

   void ASCIIRWriter::postWrite()
   {
      if(this->doesIO())
      {
         // Close the file
         this->close();
      }
   }

}
