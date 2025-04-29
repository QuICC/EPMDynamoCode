/** \file CFLTimestepFile.cpp
 *  \brief Source of the time information file
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/ASCII/CFLTimestepFile.hpp"

// Project includes
//
#include "IO/ASCII/CFLTimestepFileDefs.hpp"

namespace EPMDynamo {

   CFLTimestepFile::CFLTimestepFile(std::string name, const TimestepParameters &tsParams)
      : ASCIIEWriter(name + CFLTimestepFileDefs::BASENAME, CFLTimestepFileDefs::EXTENSION, CFLTimestepFileDefs::HEADER, CFLTimestepFileDefs::VERSION), mrTSParams(tsParams)
   {
   }

   void CFLTimestepFile::write()
   {
      if(this->doesIO())
      {
         //Do pre write operations
         this->preWrite();

         // Write timestep
         this->mFile << this->mrTSParams.time() << "  " << this->mrTSParams.dt() << "   " << this->mrTSParams.globalCFLs().transpose() << "   " << this->mrTSParams.rtpCFLs().transpose() <<  "   " << this->mrTSParams.errCFL() << "   " << this->mrTSParams.rtpCFLPos().transpose() <<  std::endl;

         //Do pre write operations
         this->postWrite();
      }
   }
}
