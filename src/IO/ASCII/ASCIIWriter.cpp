/** \file ASCIIWriter.cpp
 *  \brief Source of the general ASCII writer
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/ASCII/ASCIIWriter.hpp"

// Project includes
//
#include "General/EPMxception.hpp"

namespace EPMDynamo {

   ASCIIWriter::ASCIIWriter(std::string name, std::string ext, std::string header, std::string version)
      : ASCIIFile(name, ext, header, version)
   {
   }

   void ASCIIWriter::open()
   {
      if(this->doesIO())
      {
         this->mFile.open(this->filename().c_str());

         if(! this->mFile.is_open())
         {
            throw EPMException("ASCIIWriter::open", "Couldn't open file " + this->filename() + "!");
         }
      }
   }

   void ASCIIWriter::close()
   {
      if(this->doesIO())
      {
         this->mFile.close();
      }
   }

}
