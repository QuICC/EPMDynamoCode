/** \file SplittingFile.cpp
 *  \brief Source of the splitting information file
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/ASCII/SplittingFile.hpp"

// Project includes
//
#include "IO/ASCII/SplittingFileDefs.hpp"

namespace EPMDynamo {

   SplittingFile::SplittingFile(const SmartTruncation pTrunc, std::string name)
      : ASCIIRWriter(name + SplittingFileDefs::BASENAME, SplittingFileDefs::EXTENSION, SplittingFileDefs::HEADER, SplittingFileDefs::VERSION), mpTrunc(pTrunc)
   {
   }

}
