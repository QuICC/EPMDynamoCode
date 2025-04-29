/** \file FDSHSplittingFile.cpp
 *  \brief Source of the FDSH space splitting information file
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/ASCII/FDSHSplittingFile.hpp"

// Project includes
//

namespace EPMDynamo {

   FDSHSplittingFile::FDSHSplittingFile(const SmartTruncation pTrunc)
      : SplittingFile(pTrunc, "FDSH")
   {
   }

   void FDSHSplittingFile::write()
   {
      if(this->doesIO())
      {
         // Do pre write operations
         this->preWrite();

         // Write header information
         this->mFile << "# Nr=" << this->mpTrunc->sim()->rad()->nR() << std::endl;
         this->mFile << "# L=" << this->mpTrunc->sim()->hoz()->maxL() << std::endl;
         this->mFile << "# M=" << this->mpTrunc->sim()->hoz()->maxM() << std::endl;
         this->mFile << "Nr" << "\t\t" << "M" << "\t\t" << "cpuID" << std::endl;

         // Loop over all CPUs
         for(int id=0; id < this->mpTrunc->para().nCore();++id)
         {
            for(int iM=0; iM < this->mpTrunc->remote(id)->fdsh()->nM();++iM)
            {
               for(int n=0;n < this->mpTrunc->remote(id)->fdsh()->nR(iM); ++n)
               {
                     this->mFile << this->mpTrunc->remote(id)->fdsh()->r0(iM)+n << "\t\t";
                     this->mFile << this->mpTrunc->remote(id)->fdsh()->mArray()(iM) << "\t\t" << id << std::endl;
               }
            }
         }

         // Do post write operations
         this->postWrite();
      }
   }

}
