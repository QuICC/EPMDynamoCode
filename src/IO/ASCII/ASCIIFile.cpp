/** \file ASCIIFile.cpp
 *  \brief Source of the general ASCII file implementation
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/ASCII/ASCIIFile.hpp"

// Project includes
//

namespace EPMDynamo {

   ASCIIFile::ASCIIFile(std::string name, std::string ext, std::string header, std::string version)
      : FileBase(name, ext), mHeader(header), mVersion(version)
   {
   }

   const std::string ASCIIFile::FILE_HEADER = "#";

   const std::string ASCIIFile::FILE_VERSION = "#Version";
}
