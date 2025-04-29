/** \file ASCIINWriter.cpp
 *  \brief Source of the numbered ASCII file writer
 */

// System includes
//
#include <sstream>
#include <iomanip>

// External includes
//

// Class include
//
#include "IO/ASCII/ASCIINWriter.hpp"

// Project includes
//

namespace EPMDynamo {

   ASCIINWriter::ASCIINWriter(std::string name, std::string ext, std::string header, std::string version)
      : ASCIIWriter(name, ext, header, version), mCounter(0), mBaseName(name)
   {
   }

   const int ASCIINWriter::msIDWidth = 4;

   void ASCIINWriter::preWrite()
   {
      if(this->doesIO())
      {
         // Update name with counter value
         this->updateName();

         // Create file
         this->open();

         // Add header information
         this->mFile << FILE_HEADER + this->mHeader << std::endl;
         this->mFile << FILE_VERSION + " " + this->mVersion << std::endl;
      }
   }

   void ASCIINWriter::postWrite()
   {
      if(this->doesIO())
      {
         // Close the file
         this->close();

         // Increment the file counter
         ++this->mCounter;
      }
   }

   void ASCIINWriter::init()
   {
   }

   void ASCIINWriter::finalise()
   {
   }

   void ASCIINWriter::updateName()
   {
      std::ostringstream   oss;

      oss << std::setfill('0') << std::setw(msIDWidth) << this->mCounter;

      this->resetName(this->mBaseName + oss.str());
   }


}
