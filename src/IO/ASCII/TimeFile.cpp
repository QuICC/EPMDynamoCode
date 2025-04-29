/** \file TimeFile.cpp
 *  \brief Source of the time information file
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/ASCII/TimeFile.hpp"

// Project includes
//
#include "IO/ASCII/TimeFileDefs.hpp"

namespace EPMDynamo {

   TimeFile::TimeFile(std::string name, const TimestepParameters &tsParams)
      : ASCIIEWriter(name + TimeFileDefs::BASENAME, TimeFileDefs::EXTENSION, TimeFileDefs::HEADER, TimeFileDefs::VERSION), mrTSParams(tsParams)
   {
   }

   void TimeFile::write()
   {
      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         // Write timestep
         this->mFile << this->mrTSParams.time() << "  " << this->mrTSParams.dt() << std::endl;

         //Do pre write operations
         this->postWrite();
      }
   }
}
