/** \file ASCIIReader.cpp
 *  \brief Source of the general ASCII reader
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/ASCII/ASCIIReader.hpp"

// Project includes
//
#include "General/EPMException.hpp"

namespace EPMDynamo {

   ASCIIReader::ASCIIReader(std::string name, std::string ext, std::string header, std::string version)
      : ASCIIFile(name, ext, header, version)
   {
   }

   void ASCIIReader::init()
   {
      if(this->doesIO())
      {
         this->open();

         this->checkCompatibility();
      }
   }

   void ASCIIReader::open()
   {
      if(this->doesIO())
      {
         this->mFile.open(this->filename().c_str());

         if(! this->mFile.is_open())
         {
            throw EPMException("ASCIIReader::open", "Couldn't open file " + this->filename() + "!");
         }
      }
   }

   void ASCIIReader::checkCompatibility()
   {
      if(this->doesIO())
      {
         std::string fileHead;
         std::string fileVers;

         // Read header if present
         getline(this->mFile, fileHead);

         if(this->mFile.good())
         {
            // Check header
            if(fileHead.compare(FILE_HEADER + this->mHeader) == 0)
            {
               // Read version if present
               getline(this->mFile, fileVers);

               if(this->mFile.good())
               {
                  // Check version
                  if(fileVers.compare(FILE_VERSION + " " + this->mVersion) == 0)
                  {
                  } else
                  {
                     throw EPMException("ASCIIReader::checkCompatibility", "Wrong file version!");
                  }
               } else
               {
                  throw EPMException("ASCIIReader::checkCompatibility", "Missing file version!");
               }
            } else
            {
               throw EPMException("ASCIIReader::checkCompatibility", "Wrong file header!");
            }
         } else
         {
            throw EPMException("ASCIIReader::checkCompatibility", "Missing file header!");
         }
      }
   }

   void ASCIIReader::finalise()
   {
      if(this->doesIO())
      {
         this->close();
      }
   }

   void ASCIIReader::close()
   {
      if(this->doesIO())
      {
         this->mFile.close();
      }
   }

}
