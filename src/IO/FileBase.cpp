/** \file FileBase.cpp
 *  \brief Source of the basic aspect of a file
 */

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//

// External includes
//

// Class include
//
#include "IO/FileBase.hpp"

// Project includes
//

namespace EPMDynamo {

   FileBase::FileBase(std::string name, std::string ext)
      : mName(name), mExt(ext), mcIORank(0), mRank(EPMDYNAMO_RANK), mIsParallel(EPMDYNAMO_IS_PARALLEL)
   {
      if(this->mRank == this->mcIORank)
      {
         this->mDoesIO = true;
      } else
      {
         this->mDoesIO = false;
      }
   }

   void FileBase::resetName(std::string name)
   {
      this->mName = name;
   }
}
