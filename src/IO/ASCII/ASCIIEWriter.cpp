/** \file ASCIIEWriter.cpp
 *  \brief Source of the general ASCII extending writer
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "IO/ASCII/ASCIIEWriter.hpp"

// Project includes
//

namespace EPMDynamo {

   ASCIIEWriter::ASCIIEWriter(std::string name, std::string ext, std::string header, std::string version)
      : ASCIIWriter(name, ext, header, version)
   {
   }

   void ASCIIEWriter::init()
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

   void ASCIIEWriter::finalise()
   {
      if(this->doesIO())
      {
         // Close the file
         s->close();
      }
   }

   void ASCIIEWriter::preWrite()
   {
   }

   void ASCIIEWriter::postWrite()
   {
   }

}
